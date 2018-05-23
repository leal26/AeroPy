from scipy.interpolate import interp1d
import numpy as np
import math

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
    
def CST_3D(Bu, Bl, span, N={'eta':[0,1], 'N1':[.5, .5], 'N2':[1., 1.], 'chord':[1., 0]},
           mesh = (100,100), chord = {'eta':[0,1], 'A':[1.], 'N1':1, 'N2':1, 'initial':1., 'final':0.1}, 
           sweep = {'eta':[0,1], 'A':[1.], 'N1':1, 'N2':1, 'initial':0, 'final':0},
           twist = {'eta':[0,1], 'A':[1.], 'N1':1, 'N2':1, 'initial':0, 'final':1},
           surfaces = 'both'):
    """
    - Bu: upper shape coefficients
    - Bl: lower shape coefficients
    - mesh: list of number of points in x and y
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

        Nx = len(B)-1
        Ny = len(B[0])-1

        output = 0
        for i in range(Nx+1):
            for j in range(Ny+1):
                output += B[i][j]*S_i(i, Nx, psi)*S_i(j, Ny, eta)
        return output

    def C(psi, eta):
        """Class function"""
        psi_max = N1(eta)/(N1(eta)+N2(eta))
        C_max = 2*((psi_max)**N1(eta))*((1.-psi_max)**N2(eta))
        output = ((psi)**N1(eta))*((1.-psi)**N2(eta))/C_max/2
        return output

    def Rotation(ux, uy, uz, theta):
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
    # Define non-dimensional domains

    psi = np.linspace(0,1,mesh[1])
    eta = np.linspace(0,1,mesh[1])
    zeta_u = np.zeros(mesh)
    zeta_l = np.zeros(mesh)
    
    # Interpolate class function coefficients

    N1 = interp1d(N['eta'], N['N1'])
    N2 = interp1d(N['eta'], N['N2'])
    

    for i in range(mesh[0]):
        for j in range(mesh[1]):
            if surfaces == 'upper' or surfaces == 'both':
                zeta_u[j][i] = C(psi[i], eta[j])*S(Bu, psi[i], eta[j])
            if surfaces == 'lower' or surfaces == 'both':
                zeta_l[j][i] = -C(psi[i], eta[j])*S(Bl, psi[i], eta[j])

    chord_distribution = CST(eta, chord['eta'][1], chord['initial'], Au=chord['A'], N1=chord['N1'], N2=chord['N2'], deltasLE=chord['final'])
    sweep_distribution = CST(eta, sweep['eta'][1], deltasz = sweep['final'], Au=sweep['A'], N1=sweep['N1'], N2=sweep['N2'])
    chord_distribution = chord_distribution[::-1]
    sweep_distribution = sweep_distribution
    twist_distribution = CST(eta, twist['eta'][1], twist['initial'], Au=twist['A'], N1=twist['N1'], N2=twist['N2'], deltasLE=twist['final'])

    # taper_function(eta, shape = 'linear', N)
    x = np.zeros(len(psi))
    for i in range(len(x)):
        x[i] = psi[i]*chord_distribution[i]
    y = eta

    X = np.zeros(mesh)
    Y = np.zeros(mesh)
    X_u = np.zeros(mesh)
    X_l = np.zeros(mesh)
    Y_u = np.zeros(mesh)
    Y_l = np.zeros(mesh)
    Z_u = np.zeros(mesh)
    Z_l = np.zeros(mesh)
    for i in range(mesh[0]):
        for j in range(mesh[1]):
            X[j][i] = psi[i]*chord_distribution[j] + sweep_distribution[j]
            Y[j][i] = span*eta[j]
            Z_u[j][i] = zeta_u[j][i] *chord_distribution[j]
            Z_l[j][i] = zeta_l[j][i] *chord_distribution[j]
            if twist != None:
                xr = twist['psi_root']
                zr = twist['zeta_root']
                R = Rotation(twist['axis'][0], twist['axis'][1], twist['axis'][2],
                         twist_distribution[j])
                if surfaces == 'upper' or surfaces == 'both':
                    t_u = np.array([X[j][i]-xr,Y[j][i],Z_u[j][i]-zr])
                    T_u = R.dot(t_u)
                    X_u[j][i] = xr + T_u[0]
                    Y_u[j][i] = T_u[1]
                    Z_u[j][i] = zr + T_u[2]
                    
                if surfaces == 'lower' or surfaces == 'both':
                    t_l = np.array([X[j][i]-xr,Y[j][i],Z_l[j][i]-zr])
                    T_l = R.dot(t_l)
                    X_l[j][i] = xr + T_l[0]
                    Y_l[j][i] = T_l[1]
                    Z_l[j][i] = zr + T_l[2]
    # Formating data for output
    if surfaces == 'upper' or surfaces == 'both':
        output_u = {}
        if twist == None:
            output_u['x'] = X
            output_u['y'] = Y
        else:
            output_u['x'] = X_u
            output_u['y'] = Y_u
        output_u['z'] = Z_u

    if surfaces == 'upper' or surfaces == 'both':
        output_l = {}
        if twist == None:
            output_l['x'] = X
            output_l['y'] = Y
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
    # One of the diameters
    initial_chord = 2.
    final_chord = 1.
    # Nosecone height
    span = 4.
    # Shape coefficient for cross section (if A=1, circular, otherwise it is an ellipse)
    A = .5
    # location of the nosecone tip
    initial_nosecone_x = 0
    final_nosecone_x = 2.
    # Class coefficient for chord distribution (Nb=.5, elliptical, Nb=1, Haack series)
    Nb = 0.0
 
    axis = [final_nosecone_x - initial_nosecone_x, span,0]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    B = [[A], [A]]
    #B = [[A], [A]]
    Na = 0.0
    mesh =(20,20)
    N={'eta':[0,0.5,.8,1], 'N1':[.5, .5, 1., 1.], 'N2':[1., 1., 1., 1.]}
    chord = {'eta':[0,1], 'A':[0.], 'N1':Na, 'N2':Nb, 
             'initial':initial_chord,'final':final_chord}
    sweep = {'eta':[0,1], 'A':[0.], 'N1':Nb, 'N2':Na, 
             'initial':initial_nosecone_x, 'final':final_nosecone_x}
    twist = {'eta':[0,1], 'A':[0.], 'N1':1, 'N2':1,'initial':0.5, 'final':0, 'psi_root':0.2,'zeta_root':0, 'axis':axis}    
    
    output = CST_3D(B, B,  span, N, mesh, chord, sweep, twist)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    for surface in output:
        ax.plot_surface(output[surface]['x'], output[surface]['z'], 
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