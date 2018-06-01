from scipy.interpolate import interp1d
import numpy as np
import math
import pickle

from aeropy.geometry.airfoil import CST, create_x

def interpolation_function(eta, shape = 'linear', 
                           points = {'eta':[0,1], 'to_interpolate':[1,.7]}):
    """Calculate 'to_interpolate' along span of the wing.

    - If linear, interpolation function is a conjuction of lines connecting points
    - Possible shapes are the same as interp1d: ('linear', 'nearest',
        'zero', 'slinear', 'quadratic', 'cubic' where 'zero', 'slinear',
         'quadratic' and 'cubic' refer to a spline interpolation of zeroth,
          first, second or third order"""
    function = interp1d(points['eta'], points['to_interpolate'])
    return function(eta)
    
def CST_3D(B, span=1., mesh = (10,10), eta_control = [0, 1], N1_control = [.5, .5],
           N2_control = [1., 1.], chord_control = [1., 1.], sweep_control = [0., 0.], 
           twist_control = [0., 0.], shear_control = [0.,0.], twist_axis = {'x':0, 'z':0, 'vector':[0,1,0]},
           surfaces = 'both'):
    """
    - B: input that defines shape coefficients:
         - if surfaces == 'lower': B = Bl (upper shape coefficients)
         - if surfaces == 'upper': B = Bu (lower shape coefficients)
         - if type==list and surfaces=='both': B=[Bu, Bl]

    - mesh: input for mesh. If type(mesh) == list, a uniform mesh generator
            is used and mesh =  list of number of points in x and y. If len==2,
            the upper and lower surfaces have the same mesh. If len==4, they are
            different. To import a mesh, use a dictionary with the following
            formatting:
            mesh = {'psi_u':psi_u, 'eta_u':eta_u, 'Neta_u':Npsi_u, 
                    'Neta_u':Neta_u,
                    'psi_l':psi_l, 'eta_l':eta_l, 'Neta_l':Npsi_l, 
                    'Neta_u':Neta_l}
    - surfaces: surfaces to be analyzed (options: 'upper', 'lower', and 'both'.
    """
    def S(B, psi, eta):
        """ Cross section shape function. Validated for high dimensions.
           To debug just verify if it turns all ones when B=ones"""
        def S_i(r, n, psi):
            """Shape function"""
            value = K(r,n)*(psi**r)*(1.-psi)**(n-r)
            return value

        # Bersntein Polynomial
        def K(r,n):
            K=math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
            return K

        Ny = len(B)-1
        Nx = len(B[0])-1

        output = 0
        for i in range(Nx+1):
            for j in range(Ny+1):
                output += B[j][i]*S_i(i, Nx, psi)*S_i(j, Ny, eta)
        return output

    def C(psi, eta):
        """Class function"""
        psi_max = N1(eta)/(N1(eta)+N2(eta))
        C_max = 2*((psi_max)**N1(eta))*((1.-psi_max)**N2(eta))
        output = ((psi)**N1(eta))*((1.-psi)**N2(eta))/C_max/2
        return output

    def Rotation(euler_vector, theta):
        ux = euler_vector[0]
        uy = euler_vector[1]
        uz = euler_vector[2]
        # Check if vector are normalized
        norm = math.sqrt(ux**2+uy**2+uz*2)
        if norm != 1:
            ux = ux/norm
            uy = uy/norm
            uz = uz/norm
        R11 = math.cos(theta) + ux**2*(1-math.cos(theta))
        R12 = ux*uy*(1-math.cos(theta)) - uz*math.sin(theta)
        R13 = ux*uz*(1-math.cos(theta)) + uy*math.sin(theta)
        R21 = ux*uy*(1-math.cos(theta)) + uz*math.sin(theta)
        R22 = math.cos(theta) + uy**2*(1-math.cos(theta))
        R23 = uy*uz*(1-math.cos(theta)) - ux*math.sin(theta)
        R31 = ux*uz*(1-math.cos(theta)) - uy*math.sin(theta)
        R32 = uy*uz*(1-math.cos(theta)) + ux*math.sin(theta)
        R33 = math.cos(theta) + uz**2*(1-math.cos(theta))
        R = np.array([[R11, R12, R13],
                      [R21, R22, R23],
                      [R31, R32, R33]])
        return R

    def uniform_mesh_generator(mesh):
        '''Default genrator for uniform meshes. If 
           len(mesh)==2, upper and lower have same
           mesh'''
        # 
        if len(mesh) == 2:
            Npsi_u = mesh[0]
            Neta_u = mesh[1]
            Npsi_l = mesh[0]
            Neta_l = mesh[1]
        elif len(mesh) == 4:
            Npsi_u = mesh[0]
            Neta_u = mesh[1]
            Npsi_l = mesh[2]
            Neta_l = mesh[3]
        psi_u, eta_u = np.meshgrid(np.linspace(0,1,Npsi_u), 
                                   np.linspace(0,1,Neta_u))
        psi_l, eta_l = np.meshgrid(np.linspace(0,1,Npsi_l), 
                                   np.linspace(0,1,Neta_l))

        mesh = {'psi_u':psi_u, 'eta_u':eta_u, 'Npsi_u':Npsi_u, 
                'Neta_u':Neta_u,
                'psi_l':psi_l, 'eta_l':eta_l, 'Npsi_l':Npsi_l, 
                'Neta_l':Neta_l}
        return mesh

    # Check if inputs are formatted properly
    for property in [N1_control, N2_control, chord_control, 
                     sweep_control, twist_control, shear_control]:
        if len(eta_control) != len(property):
            raise('All geometric properties must have same length')

    # Processing inputs
    if surfaces == 'both':
        Bu = B[0]
        Bl = B[1]
    elif surfaces == 'upper':
        Bu = B
    elif surfaces == 'lower':
        Bl = B

    if type(mesh) != dict:
        mesh = uniform_mesh_generator(mesh)
        psi_u = mesh['psi_u']
        psi_l = mesh['psi_l']
        eta_u = mesh['eta_u']
        eta_l = mesh['eta_l']
        Npsi_u = mesh['Npsi_u']
        Npsi_l = mesh['Npsi_l']
        Neta_u = mesh['Neta_u']
        Neta_l = mesh['Neta_l']

    # Define non-dimensional domains
    zeta_u = np.zeros((Neta_u, Npsi_u))
    zeta_l = np.zeros((Neta_l, Npsi_l))

    # Interpolate class function coefficients and other properties
    N1 = interp1d(eta_control, N1_control)
    N2 = interp1d(eta_control, N2_control)
    chord = interp1d(eta_control, chord_control)
    sweep = interp1d(eta_control, sweep_control)
    twist = interp1d(eta_control, twist_control)   
    shear = interp1d(eta_control, shear_control)

    # Calculating non-dimensional surface
    if surfaces == 'upper' or surfaces == 'both':
        for i in range(Npsi_u):
            for j in range(Neta_u):
                zeta_u[j][i] = C(psi_u[j][i], eta_u[j][i]
                         )*S(Bu, psi_u[j][i], eta_u[j][i])

    if surfaces == 'lower' or surfaces == 'both':
        for i in range(Npsi_l):
            for j in range(Neta_l):
                zeta_l[j][i] = -C(psi_l[j][i], eta_l[j][i]
                          )*S(Bl, psi_l[j][i], eta_l[j][i])

    # Calculate surface on physical domain
    X_u = np.zeros((Neta_u, Npsi_u))
    X_l = np.zeros((Neta_l, Npsi_l))
    Y_u = np.zeros((Neta_u, Npsi_u))
    Y_l = np.zeros((Neta_l, Npsi_l))
    Z_u = np.zeros((Neta_u, Npsi_u))
    Z_l = np.zeros((Neta_l, Npsi_l))
    if surfaces == 'upper' or surfaces == 'both':
        for i in range(Npsi_u):
            for j in range(Neta_u):
                X_u[j][i] = psi_u[j][i]*chord(eta_u[j][i]) + sweep(eta_u[j][i])
                Y_u[j][i] = span*eta_u[j][i]
                Z_u[j][i] = zeta_u[j][i]*chord(eta_u[j][i])

                # Twisting component
                R = Rotation(twist_axis['vector'], twist(eta_u[j][i]))
                t_u = np.array([X_u[j][i] - twist_axis['x'],
                                Y_u[j][i],
                                Z_u[j][i] - twist_axis['z']])
                T_u = R.dot(t_u)
                X_u[j][i] = twist_axis['x'] + T_u[0]
                Y_u[j][i] = T_u[1]
                Z_u[j][i] = twist_axis['z'] + T_u[2] + shear(eta_l[j][i])

    if surfaces == 'lower' or surfaces == 'both':
        for i in range(Npsi_l):
            for j in range(Neta_l):
                X_l[j][i] = psi_l[j][i]*chord(eta_l[j][i]) + sweep(eta_l[j][i])
                Y_l[j][i] = span*eta_l[j][i]
                Z_l[j][i] = zeta_l[j][i]*chord(eta_l[j][i])

                # Twisting component
                R = Rotation(twist_axis['vector'], twist(eta_u[j][i]))
                t_l = np.array([X_l[j][i] - twist_axis['x'],
                                Y_l[j][i],
                                Z_l[j][i] - twist_axis['z']])
                T_l = R.dot(t_l)
                X_l[j][i] = twist_axis['x'] + T_l[0]
                Y_l[j][i] = T_l[1]
                Z_l[j][i] = twist_axis['z'] + T_l[2] + shear(eta_l[j][i])

    # Formating data for output
    if surfaces == 'upper' or surfaces == 'both':
        output_u = {}
        if twist == None:
            output_u['x'] = X_u
            output_u['y'] = Y_u
        else:
            output_u['x'] = X_u
            output_u['y'] = Y_u
        output_u['z'] = Z_u

    if surfaces == 'upper' or surfaces == 'both':
        output_l = {}
        if twist == None:
            output_l['x'] = X_l
            output_l['y'] = Y_l
        else:
            output_l['x'] = X_l
            output_l['y'] = Y_l
        output_l['z'] = Z_l
    if surfaces == 'upper':
        return output_u
    elif surfaces == 'lower':
        return output_l
    else:
        return {'upper':output_u,
                'lower':output_l}

if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm

    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Inputs
    span = 4.
    # Shape coefficient for cross section (if A=1, circular, otherwise it is an ellipse)
    Au = np.array([0.172802, 0.167353, 0.130747, 0.172053, 0.112797, 0.168891])
    Al = np.array([0.163339, 0.175407, 0.134176, 0.152834, 0.133240, 0.161677])

    B = [[Au,Au], [Al,Al]]
    output = CST_3D(B, chord_control = [1., .1], twist_control = [0, .1])    
    # A = 0.5
    # # Mesh size in case default mesh generator is used
    # mesh = (20,40)
    # B = [[A], [A]] 
    # output = CST_3D([B,B], span =4., shear = [0,.5], twist = [0,0.1],
                    # chord = [1.0,.1])
    pickle.dump( output, open( "test_case.p", "wb" ) )
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    for surface in output:
        ax.scatter(output[surface]['x'], output[surface]['z'], 
                        output[surface]['y'], cmap=plt.get_cmap('jet'),
                        linewidth=0, antialiased=False)
    # cset = ax.contour(X, Z_u, Y, zdir='z', offset=0, cmap=cm.coolwarm)
    # cset = ax.contour(X, Z_l, Y, zdir='z', offset=0,  cmap=cm.coolwarm)
    # cset = ax.contour(X, Z_u, Y, zdir='x', offset=-.1, cmap=cm.coolwarm)
    # cset = ax.contour(X, Z_l, Y, zdir='x', offset=-.1, cmap=cm.coolwarm)
    # cset = ax.contour(X, Z_u, Y, zdir='y', offset =0.5,  cmap=cm.coolwarm)
    # cset = ax.contour(X, Z_l, Y, zdir='y', offset =0.5,  cmap=cm.coolwarm)
    
    # Customize the z axis.
    ax.set_zlim(0, 4)

    max_range = np.array([output['upper']['x'].max()-output['upper']['x'].min(),
                          output['lower']['y'].max()-output['lower']['y'].min(),
                          output['upper']['z'].max()-output['lower']['z'].min(),
                          output['upper']['y'].max()-output['upper']['y'].min(),
                          output['lower']['y'].max()-output['lower']['y'].min()]).max() / 2.0

    mid_x = (output['upper']['x'].max()+output['lower']['x'].min()) * 0.5
    mid_y = (output['upper']['y'].max()+output['lower']['x'].min()) * 0.5
    mid_z = (output['upper']['z'].max()+output['lower']['z'].min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_z - max_range, mid_z + max_range)
    ax.set_zlim(mid_y - max_range, mid_y + max_range)
    plt.xlabel('x')
    plt.ylabel('z')

    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    
    for surface in output:
        ax.scatter(output[surface]['x'], output[surface]['z'])

    plt.xlabel('x')
    plt.ylabel('z')
    plt.show()

    # fig = plt.figure()
    # ax = fig.gca(projection='3d')

    # ax.plot_trisurf(X_u.flatten(),  Y_u.flatten(), Z_u.flatten(), linewidth=0.2, antialiased=True,)
    # ax.plot_trisurf(X_l.flatten(),  Y_l.flatten(), Z_l.flatten(), linewidth=0.2, antialiased=True,)
    # ax.plot([sweep['initial'] + twist['psi_root']*chord['initial'],sweep['initial'] + twist['psi_root']*chord['initial'] + axis[0]],[0,axis[1]])
    # ax.set_xlim(mid_x - max_range, mid_x + max_range)
    # ax.set_ylim(mid_y - max_range, mid_y + max_range)
    # ax.set_zlim(mid_z - max_range, mid_z + max_range)
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.show()