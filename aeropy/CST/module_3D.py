from scipy.interpolate import interp1d
import numpy as np
import math
import pickle

from aeropy.geometry.airfoil import CST, create_x
from aeropy.geometry.wing import Rotation_euler_vector
    
class control_points():
    def __init__(self):
        self.eta = [0,1]
        self.N1 = [.5, .5]
        self.N2 = [1., 1.]
        self.chord = [1.,1.]
        self.sweep = [0.,0.]
        self.twist = [0., 0.]
        self.shear = [0., 0.]
        self.twist_axis = [0,1,0]
        self.twist_origin = {'x':0, 'z':0}
        
        self.half_span = 1.

    def check_attributes(self):
        # Check if inputs are formatted properly
        for property in ['N1', 'N2', 'chord', 'sweep', 'twist', 'shear']:
            if len(self.eta) != len(getattr(self,property)):
                raise Exception('All geometric properties must have same length')   
    def set_functions(self):
        self.fN1 = interp1d(self.eta, self.N1)
        self.fN2 = interp1d(self.eta, self.N2)
        self.fchord = interp1d(self.eta, self.chord)
        self.fsweep = interp1d(self.eta, self.sweep)
        self.ftwist = interp1d(self.eta, self.twist)
        self.fshear = interp1d(self.eta, self.shear)

def CST_3D(B, mesh = (10,10), cp = control_points(),
           surfaces = 'both', origin = {'x':0, 'y':0, 'z':0}, axis_order = ['x','y','z']):
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
            mesh = {'psi_u':psi_u, 'eta_u':eta_u, 'Npsi_u':Npsi_u, 
                    'Neta_u':Neta_u,
                    'psi_l':psi_l, 'eta_l':eta_l, 'Npsi_l':Npsi_l, 
                    'Neta_u':Neta_l}
            Meshes can also be imported in physical space as follows
            mesh = {'x_u', 'y_u', 'x_l', 'y_l'}
            Furthermore, you can define for only lower or upper surface.
            No need to have both.
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

        psi_max = cp.fN1(eta)/(cp.fN1(eta)+cp.fN2(eta))
        C_max = 2*((psi_max)**cp.fN1(eta))*((1.-psi_max)**cp.fN2(eta))
        output = ((psi)**cp.fN1(eta))*((1.-psi)**cp.fN2(eta))/C_max/2

        return output


    
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

    def process_mesh(mesh, cp = cp):
        if type(mesh) != dict:
                mesh = uniform_mesh_generator(mesh)
        elif 'x_u' in mesh:
            input = {'x': [x - origin['x'] for x in mesh['x_u']],
                     'y': [y - origin['y'] for y in mesh['y_u']]}
            upper = nondimensionalize(input, cp)
            mesh['psi_u'] = upper['psi']
            mesh['eta_u'] = upper['eta']
            mesh['Npsi_u'] = len(upper['psi'])
            mesh['Neta_u'] = len(upper['eta'])
        if 'x_l' in mesh:   
            input = {'x': [x - origin['x'] for x in mesh['x_l']],
                     'y': [y - origin['y'] for y in mesh['y_l']]} 
            
            lower = nondimensionalize(input, cp)
            mesh['psi_l'] = lower['psi']
            mesh['eta_l'] = lower['eta']
            mesh['Npsi_l'] = len(lower['psi'])
            mesh['Neta_l'] = len(lower['eta'])
        if 'x' in mesh:
            input = {}
            for key in mesh:
                input[key] = [x - origin[key] for x in mesh[key]]

            coordinates = nondimensionalize({'x':input[axis_order[0]], 'y':input[axis_order[1]]},
                                      cp=cp)

            if surfaces == 'upper' or surfaces == 'both':
                mesh['psi_u'] = [coordinates['psi']]
                mesh['eta_u'] = [coordinates['eta']]
                mesh['Npsi_u'] = len(coordinates['psi'])
                mesh['Neta_u'] = 1
            if surfaces == 'lower' or surfaces == 'both':
                mesh['psi_l'] = [coordinates['psi']]
                mesh['eta_l'] = [coordinates['eta']]
                mesh['Npsi_l'] = len(coordinates['psi'])
                mesh['Neta_l'] = 1
        if 'psi' in mesh:
            if 'y' in mesh:
                mesh['eta'] = [y/cp.half_span for y in mesh['y']]
            if surfaces == 'upper' or surfaces == 'both':
                mesh['psi_u'] = [mesh['psi']]
                mesh['eta_u'] = [mesh['eta']]
                mesh['Npsi_u'] = len(mesh['psi'])
                mesh['Neta_u'] = 1
            if surfaces == 'lower' or surfaces == 'both':
                mesh['psi_l'] = [mesh['psi']]
                mesh['eta_l'] = [mesh['eta']]
                mesh['Npsi_l'] = len(mesh['psi'])
                mesh['Neta_l'] = 1               
        return mesh
        
    cp.check_attributes()
    cp.set_functions()

    # Processing inputs
    if surfaces == 'both':
        Bu = B[0]
        Bl = B[1]
    elif surfaces == 'upper':
        Bu = B
    elif surfaces == 'lower':
        Bl = B

    # if z coordinate defined do inverse problem instead of calculating
    # everything again
    if 'zeta' in mesh:
        psi_list = []
        eta_list = []
        zeta_list = []
        # If zeta is a list it will work. Otherwise just say length 1
        try:
            Nzeta = len(mesh['zeta'])
        except:
            Nzeta = 1
        for i in range(Nzeta):
            input = {}
            for key in mesh:
                try:
                    input[key] = mesh[key][i]
                except:
                    input[key] = mesh[key]
            psi, eta, zeta = inverse_problem(coordinates = input, 
                                 B = B, cp=cp, 
                                 surface = surfaces)
            psi_list.append(psi)
            eta_list.append(eta)
            zeta_list.append(zeta)
        if surfaces == 'upper' or surfaces == 'both':
            psi_u = [psi_list]
            eta_u = [eta_list]
            Npsi_u = 1
            Neta_u = 1
            zeta_u = [zeta_list]
        if surfaces == 'lower' or surfaces == 'both':
            psi_l = [psi_list]
            eta_l = [eta_list]
            Npsi_l = 1
            Neta_l = 1
            zeta_l = [zeta_list]  
    else:
        # Process mesh
        mesh = process_mesh(mesh, cp)
                        
        # Calculating non-dimensional surface
        if surfaces == 'upper' or surfaces == 'both':
            # define non-dimensional domains
            psi_u = mesh['psi_u']
            eta_u = mesh['eta_u']
            Npsi_u = mesh['Npsi_u']
            Neta_u = mesh['Neta_u']
            zeta_u = np.zeros((Neta_u, Npsi_u))
            for i in range(Npsi_u):
                for j in range(Neta_u):
                    zeta_u[j][i] = C(psi_u[j][i], eta_u[j][i]
                             )*S(Bu, psi_u[j][i], eta_u[j][i])

        if surfaces == 'lower' or surfaces == 'both':
            # define non-dimensional domains
            psi_l = mesh['psi_l']
            eta_l = mesh['eta_l']
            Npsi_l = mesh['Npsi_l']
            Neta_l = mesh['Neta_l']
            zeta_l = np.zeros((Neta_l, Npsi_l))
            for i in range(Npsi_l):
                for j in range(Neta_l):
                    zeta_l[j][i] = -C(psi_l[j][i], eta_l[j][i]
                              )*S(Bl, psi_l[j][i], eta_l[j][i])

    # Calculate surface on physical domain
    if surfaces == 'upper' or surfaces == 'both':
        X_u = np.zeros((Neta_u, Npsi_u))
        Y_u = np.zeros((Neta_u, Npsi_u))
        Z_u = np.zeros((Neta_u, Npsi_u))
        for i in range(Npsi_u):
            for j in range(Neta_u):
                X_u[j][i] = psi_u[j][i]*cp.fchord(eta_u[j][i]) + cp.fsweep(eta_u[j][i])
                Y_u[j][i] = cp.half_span*eta_u[j][i]
                Z_u[j][i] = zeta_u[j][i]*cp.fchord(eta_u[j][i])

                # Twisting component
                R = Rotation_euler_vector(cp.twist_axis, cp.ftwist(eta_u[j][i]))
                t_u = np.array([X_u[j][i] - cp.twist_origin['x'],
                                Y_u[j][i],
                                Z_u[j][i] - cp.twist_origin['z']])
                T_u = R.dot(t_u)
                X_u[j][i] = cp.twist_origin['x'] + T_u[0]
                Y_u[j][i] = T_u[1]
                Z_u[j][i] = cp.twist_origin['z'] + T_u[2] + cp.fshear(eta_u[j][i])

    if surfaces == 'lower' or surfaces == 'both':
        X_l = np.zeros((Neta_l, Npsi_l))
        Y_l = np.zeros((Neta_l, Npsi_l))
        Z_l = np.zeros((Neta_l, Npsi_l))
        for i in range(Npsi_l):
            for j in range(Neta_l):
                X_l[j][i] = psi_l[j][i]*cp.fchord(eta_l[j][i]) + cp.fsweep(eta_l[j][i])
                Y_l[j][i] = cp.half_span*eta_l[j][i]
                Z_l[j][i] = zeta_l[j][i]*cp.fchord(eta_l[j][i])

                # Twisting component
                R = Rotation_euler_vector(cp.twist_axis, cp.ftwist(eta_l[j][i]))
                t_l = np.array([X_l[j][i] - cp.twist_origin['x'],
                                Y_l[j][i],
                                Z_l[j][i] - cp.twist_origin['z']])
                T_l = R.dot(t_l)
                X_l[j][i] = cp.twist_origin['x'] + T_l[0]
                Y_l[j][i] = T_l[1]
                Z_l[j][i] = cp.twist_origin['z'] + T_l[2] + cp.fshear(eta_l[j][i])

    # Formating data for output
    if surfaces == 'upper' or surfaces == 'both':
        output_u = {}
        output_u[axis_order[0]] = X_u
        output_u[axis_order[1]] = Y_u
        output_u[axis_order[2]] = Z_u
        for i in range(len(axis_order)):
            output_u[axis_order[i]] += origin[axis_order[i]]
            
    if surfaces == 'lower' or surfaces == 'both':
        output_l = {}
        output_l[axis_order[0]] = X_l
        output_l[axis_order[1]] = Y_l
        output_l[axis_order[2]] = Z_l
        for i in range(len(axis_order)):
            output_l[axis_order[i]] += origin[axis_order[i]]
    if surfaces == 'upper':
        return output_u
    elif surfaces == 'lower':
        return output_l
    else:
        return {'upper':output_u,
                'lower':output_l}

def nondimensionalize(Raw_data, cp = control_points()):

    cp.check_attributes()
    cp.set_functions()
    
    # Defining coordinates
    x = Raw_data['x']
    y = Raw_data['y']
    if 'z' in Raw_data:
        z = Raw_data['z']
    else:
        z = np.zeros(len(x))

    # Converting to x0-y0-z0 domain
    x0 = []
    y0 = []
    z0 = []

    eta = [a/cp.half_span for a in y]

    twist = cp.ftwist(eta)
    shear = cp.fshear(eta)
    
    for i in range(len(y)):
        # transpose/inverse of rotation matrix
        R_i = np.transpose(Rotation_euler_vector(cp.twist_axis, twist[i])) 
        u = np.array([x[i]-cp.twist_origin['x'],
                      y[i],
                      z[i]-shear[i]])
        u0 = R_i.dot(u)
        
        x0.append(u0[0] + cp.twist_origin['x'])
        y0.append(u0[1])
        z0.append(u0[2])
    # Converting to psi-eta-xi domain
    [psi, eta, xi] = uncoupled_domain(x0, z0, y0, cp=cp,  inverse = True)
    return {'psi': psi, 'eta': eta, 'xi': xi}

def uncoupled_domain(psi, xi, eta, cp = control_points(), inverse = False):
    '''Converts from Parameterized domain to uncoupled domain (not rotated).
       - psi, xi, eta are the coordintes of the points to be analyzed
       - eta_control points where we know values such as chord etc'''

    cp.check_attributes()
    cp.set_functions()
    
    if inverse:
        eta = [a/cp.half_span for a in eta]
    chord = cp.fchord(eta)
    sweep = cp.fsweep(eta)
    x0 = []
    y0 = []
    z0 = []
    for i in range(len(eta)):
        if not inverse:
            x0.append(chord[i]*psi[i] + sweep[i])
            y0.append(cp.half_span*eta[i])
            z0.append(chord[i]*xi[i])
        else:
            x0.append((psi[i] - sweep[i])/chord[i])
            y0.append(eta[i])
            z0.append(xi[i]/chord[i])
    return x0,y0,z0    

def inverse_problem(coordinates = {'psi':0.5, 'zeta':0.5}, B = [[1]], 
                    cp=control_points(), surface = 'upper'):
    from scipy.optimize import fixed_point, minimize
    
    def S_modified(B, psi, eta):
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
                if i!=Nx or j!=Ny:
                    output += B[j][i]*S_i(i, Nx, psi)*S_i(j, Ny, eta)
                    # print B[j][i],S_i(i, Nx, psi),S_i(j, Ny, eta)s
        return output

    def C(psi, eta):
        """Class function"""

        psi_max = cp.fN1(eta)/(cp.fN1(eta)+cp.fN2(eta))
        C_max = 2*((psi_max)**cp.fN1(eta))*((1.-psi_max)**cp.fN2(eta))
        output = ((psi)**cp.fN1(eta))*((1.-psi)**cp.fN2(eta))/C_max/2

        return output
    
    
    def transcendental_equation(psi, eta, zeta, known_coordinate, surface):
        c = (zeta/C(psi, eta) - S_modified(B,psi,eta))/B[-1][-1]
        psi0 = psi
        # print 'psi', psi
        # print 'c', c
        # print 'zeta', zeta
        # print 'class', C(psi, eta)
        # print 'shape', S_modified(B,psi,eta)
        if known_coordinate == 'eta':
            psi = (c/eta**m)**(1./n)
            return psi
        elif known_coordinate == 'psi':
            eta = (c/psi**n)**(1./m)
            return eta

    def residual(psi, eta, zeta, surface):
        output = CST_3D(B, mesh = {'psi':psi, 'eta':eta}, cp = control_points(), surfaces = surface)
        return abs(zeta - output['z'])
    cp.check_attributes()
    cp.set_functions()

    m = len(B)-1
    n = len(B[0])-1
        
    if 'psi' in coordinates:
        psi = coordinates['psi']
        zeta = coordinates['zeta']
        f = lambda x: residual([psi],x,zeta, surface)
        solution = minimize(f,1., bounds = [[0,1]])
        eta = solution.x[0]
    elif 'eta' in coordinates:
        eta = coordinates['eta']
        zeta = coordinates['zeta']
        f = lambda x: residual(x,[eta],zeta,  surface)
        solution = minimize(f,1., bounds = [[0,1]])
        psi = solution.x[0]
    output = CST_3D(B, mesh = {'psi':[psi], 'eta':[eta]}, cp = control_points(), surfaces = surface)
    # Check if found actual solution
    tol = 1e-4
    if abs(output['z']-zeta) < tol:
        return psi, eta, zeta
    else:
        raise Exception('There is no point on surface with required z coordinate')
    
if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Inputs

    # Shape coefficient for cross section (if A=1, circular, otherwise it is an ellipse)
    Au = np.array([0.172802, 0.167353, 0.130747, 0.172053, 0.112797, 0.168891])
    Al = np.array([0.163339, 0.175407, 0.134176, 0.152834, 0.133240, 0.161677])
    # Wing properties (control points)
    cp = control_points()
    cp.check_attributes()
    cp.set_functions()
    B = [[Au,Au], [Al,Al]]
    output = CST_3D(B, cp = cp) 

    eta = 0.5
    psi_list = []
    eta_list = []
    zeta_list = []
    # Inverse problem
    for zeta in np.linspace(0,0.03, 10):
        [psi, eta, zeta] = inverse_problem(coordinates = {'eta':eta, 'zeta':zeta}, B = [Au,Au], 
                             cp=control_points(), surface = 'upper')

        psi_list.append(psi)
        eta_list .append(eta)
        zeta_list.append(zeta)
    # A = 0.5
    # # Mesh size in case default mesh generator is used
    # mesh = (20,40)
    # B = [[A], [A]] 
    # output = CST_3D([B,B], span =4., shear = [0,.5], twist = [0,0.1],
                    # chord = [1.0,.1])
    pickle.dump( output, open( "test_case.p", "wb" ) )
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(psi_list,  zeta_list, eta_list)
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