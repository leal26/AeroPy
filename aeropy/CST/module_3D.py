from scipy.interpolate import interp1d
import numpy as np
import math
import pickle

from aeropy.geometry.airfoil import CST, create_x
from aeropy.geometry.wing import Rotation_euler_vector
from panairwrapper.filehandling import generate_vtk_input

# class aircraft():
    # def __init__():
        # self.fuselage = 
        # self.wing = 
        # self.tail = None
    
    # def intersection():
        # return 0
    # def meshing():
        # return 0
    # def run_panair():
        # return 0
    # def run_sboom():
        # return 0
    # def run_pyldb():
        # return 0
    # def generate_vtk()
        # return 0
    # def generate_stl()
        # return 0
 
# class CST_object(surfaces, origin, axis_order):
    # def __init__():
        # self.surfaces = surfaces
        
    # def 

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

        # Prepare interpolation functions for default value
        self._set_functions()
        
    def set(self, **inputs):
        """Set control point values
        """
        for property, value in inputs.items():
            try:
                setattr(self,property,value)
            except:
                raise Exception(property + " keyword argument not recognized")

        # update values
        self._check_attributes()
        self._set_functions()

    def _check_attributes(self):
        # Check if inputs are formatted properly
        for property in ['N1', 'N2', 'chord', 'sweep', 'twist', 'shear']:
            try: 
                if len(self.eta) != len(getattr(self,property)):
                    raise Exception('All geometric properties must have same length')
            except:
                raise Exception('Excepted list. Got ' + 
                                str(type(getattr(self,property))) + 
                                ' for ' + property)
                
    def _set_functions(self):
        self.fN1 = interp1d(self.eta, self.N1)
        self.fN2 = interp1d(self.eta, self.N2)
        self.fchord = interp1d(self.eta, self.chord)
        self.fsweep = interp1d(self.eta, self.sweep)
        self.ftwist = interp1d(self.eta, self.twist)
        self.fshear = interp1d(self.eta, self.shear)
        
def CST_3D(B, mesh = (10,10), cp = control_points(),
           origin = [0,0,0], axis_order = [0,1,2],
           format = 'array', **optional_arguments):
    """
    - B: dict (keys 'upper' and 'lower') input that defines shape coefficients

    - mesh: input for mesh. If type(mesh) == list, a uniform mesh generator
            is used and mesh =  list of number of points in x and y. If len==2,
            the upper and lower surfaces have the same mesh. If len==4,      they are
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
    - 
    
    - format: controls format of output as either 'array' (size = [N,3])
              or 'matrix' (size = [N,M,3])
    
    
    """
    # Process mesh
    mesh, dimensions = process_mesh(mesh, cp, optional_arguments)
    output = {}
    for surface in B:
        # Calculating non-dimensional surface
        data_parm = parameterized_transformation(mesh[surface], 
                                 B=B[surface], surface=surface)
            
        # Dimensionalize
        data_phys= dimensionalize(data_parm, cp)
        
        # Position part in assembly_transformation
        data_assembly = assembly_transformation(data_phys, axis_order, origin)
        
        if format == 'matrix':
            data_assembly.resize(dimensions[surface] + [3])
        output[surface] = data_assembly

    if len(output.keys()) == 1:
        return output[output.keys()[0]]
    else:
        return output

def parameterized_transformation(mesh, B, surface = 'upper'):
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

        output = np.zeros(psi.shape)
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

    psi = mesh[:,0]
    eta = mesh[:,1]
    zeta = C(psi, eta)*S(B, psi, eta)
    zeta.resize([len(zeta),1])
    if surface == 'lower':
        zeta *= -1
    output = np.concatenate([mesh,zeta],1)

    return(output)
        
def dimensionalize(Raw_data, cp = control_points(), inverse = False):
    '''Converts from parameterized domain to physical domain'''
    if inverse:
        data = physical_transformation(Raw_data, cp, inverse)
        data = intermediate_transformation(data, cp, inverse)
    else:
        data = intermediate_transformation(Raw_data, cp, inverse)  
        data = physical_transformation(data, cp, inverse)
        
    return data
    

def intermediate_transformation(data_nondimensional, cp = control_points(), 
                                inverse = False):
    '''Converts from Parameterized domain to intermediate domain (not rotated).
       - psi, xi, eta are the coordinates of the points to be analyzed
       - eta_control points where we know values such as chord etc'''

    psi, eta, zeta = data_nondimensional.T
    print(psi)
    if inverse:
        eta = eta/cp.half_span
    chord = cp.fchord(eta)
    sweep = cp.fsweep(eta)

    if not inverse:
        x0 = chord*psi + sweep
        y0 = cp.half_span*eta
        z0 = chord*zeta
    else:
        x0 = (psi - sweep)/chord
        y0 = eta
        z0 = zeta/chord
    # Converting to columns
    x0.resize([len(x0),1])
    y0.resize([len(y0),1])
    z0.resize([len(z0),1])
    return np.concatenate([x0,y0,z0],1)    

def physical_transformation(data, cp = control_points(), inverse = False):
            
    if inverse:
        # Defining coordinates
        x = Raw_data['x']
        y = Raw_data['y']
        if 'z' in Raw_data:
            z = Raw_data['z']
        else:
            z = np.zeros(len(x))

        # Converting to x0-y0-z0 domain
        output = np.zeros(len(data), len(data[0]),3)

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
            
            output[i].append(u0[0] + cp.twist_origin['x'])
            y0.append(u0[1])
            z0.append(u0[2])
        # Converting to psi-eta-xi domain
        return x0, y0, z0
    else:
        eta= data[:,1]
        #displacing to center of rotation
        data[:,0] -= cp.twist_origin['x']
        data[:,1] -= cp.twist_origin['z']

        # Twisting
        R = Rotation_euler_vector(cp.twist_axis, cp.ftwist(eta))
        data = np.array([np.dot(R[:,:,l],data[l]) for l in range(len(eta))])

        #Return to position and shear
        data[:,0] -= cp.twist_origin['x']
        data[:,1] += -cp.twist_origin['z'] + cp.fshear(eta)
    return(data)

def assembly_transformation(data, axis_order, origin):
    '''Transform from physical domain to assembly domain'''
    data[:,axis_order] = data
    data += np.array(origin)
            
    return(data)

###########################
# Extra functionalities        
def inverse_problem(coordinates = {'psi':0.5, 'zeta':0.5}, B = {'upper':[[1]]}, 
                    cp=control_points()):
    from scipy.optimize import minimize
    
    def residual(psi, eta, zeta, surface):
        mesh = np.array([[psi,eta],])
        output = CST_3D(B, mesh = mesh, cp = control_points(),mesh_options=)
        return abs(zeta - output['z'])
    
    zeta = coordinates['zeta']
    if 'psi' in coordinates:
        psi = coordinates['psi']
        f = lambda x: residual([psi],x,zeta, surface)
        solution = minimize(f,1., bounds = [[0,1]])
        eta = solution.x[0]
    elif 'eta' in coordinates:
        eta = coordinates['eta']
        f = lambda x: residual(x,[eta],zeta,  surface)
        solution = minimize(f,1., bounds = [[0,1]])
        psi = solution.x[0]

    # Check if found actual solution
    tol = 1e-4
    print(solution.f)
    if abs(solution.f) < tol:
        return psi, eta, zeta
    else:
        raise Exception('There is no point on surface with required z coordinate')

###########################
# Meshing
def uniform_mesh_generator(mesh):
    '''Default genrator for uniform meshes. If 
       len(mesh)==2, upper and lower have same
       mesh'''
    # Defining dimensions matrix
    if len(mesh) == 2:
        dimensions = {'upper': list(mesh),
                      'lower': list(mesh)}
    elif len(mesh) == 4:
        dimensions = {'upper': list(mesh[:2]),
                      'lower': list(mesh[2:])}
    # Calculating grid
    upper_x, upper_y = np.meshgrid(np.linspace(0,1,dimensions['upper'][0]), 
                                   np.linspace(0,1,dimensions['upper'][1]))
    lower_x, lower_y = np.meshgrid(np.linspace(0,1,dimensions['lower'][0]), 
                                   np.linspace(0,1,dimensions['lower'][1]))

    mesh = {'upper': np.concatenate([upper_x.reshape(1,np.size(upper_x)),
                                     upper_y.reshape(1,np.size(upper_y))]).T,
            'lower': np.concatenate([lower_x.reshape(1,np.size(lower_x)),
                                     lower_y.reshape(1,np.size(lower_y))]).T}

    return mesh, dimensions

def process_mesh(mesh, cp = control_points(), mesh_options={}):
    ''' 'Universal' tool processing for following cases:
        - mesh_options not defined: mesh=[Npsi_upper, Neta_upper, Npsi_lower, Neta_lower]
                                    or mesh=[Npsi_upper, Neta_upper]=[Npsi_lower, Neta_lower]
                                    '''

    if len(mesh_options.keys())==0:
        mesh, dimensions = uniform_mesh_generator(mesh)
    elif 'x_u' in mesh:
        input = {'x': [x - origin['x'] for x in mesh['x_u']],
                 'y': [y - origin['y'] for y in mesh['y_u']]}
        upper = physical_transformation(input, cp, inverse = True)
        mesh['psi_u'] = upper['psi']
        mesh['eta_u'] = upper['eta']
        mesh['Npsi_u'] = len(upper['psi'])
        mesh['Neta_u'] = len(upper['eta'])
    if 'x_l' in mesh:   
        input = {'x': [x - origin['x'] for x in mesh['x_l']],
                 'y': [y - origin['y'] for y in mesh['y_l']]} 
        
        lower = physical_transformation(input, cp, inverse = True)
        mesh['psi_l'] = lower['psi']
        mesh['eta_l'] = lower['eta']
        mesh['Npsi_l'] = len(lower['psi'])
        mesh['Neta_l'] = len(lower['eta'])
    if 'x' in mesh:
        input = {}
        for key in mesh:
            input[key] = [x - origin[key] for x in mesh[key]]

        coordinates = physical_transformation({'x':input[axis_order[0]], 'y':input[axis_order[1]]},
                                  cp=cp, inverse=True)

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
    return mesh, dimensions
    
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
    cp.set(chord = [1.,.2],
           twist = [0, .1])
    
    B = {'upper':[Au,Au], 'lower':[Al,Al]}
    output = CST_3D(B, cp = cp, format='matrix') 

    # Generate vtk file
    generate_vtk_input(output, filename='test')
    
    eta = 0.5
    psi_list = []
    eta_list = []
    zeta_list = []
    # Inverse problem
    for zeta in np.linspace(0,0.03, 10):
        [psi, eta, zeta] = inverse_problem(coordinates = {'eta':eta, 'zeta':zeta}, B = [Au,Au], 
                             cp=control_points())

        psi_list.append(psi)
        eta_list.append(eta)
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
    # ax.scatter(psi_list,  zeta_list, eta_list)
    for surface in output:
        x,y,z = output[surface].T
        ax.scatter(x, z, y, cmap=plt.get_cmap('jet'),
                        linewidth=0, antialiased=False)
    
    # Customize the z axis.
    ax.set_zlim(0, 4)

    # max_range = np.array([output['upper']['x'].max()-output['upper']['x'].min(),
                          # output['lower']['y'].max()-output['lower']['y'].min(),
                          # output['upper']['z'].max()-output['lower']['z'].min(),
                          # output['upper']['y'].max()-output['upper']['y'].min(),
                          # output['lower']['y'].max()-output['lower']['y'].min()]).max() / 2.0

    # mid_x = (output['upper']['x'].max()+output['lower']['x'].min()) * 0.5
    # mid_y = (output['upper']['y'].max()+output['lower']['x'].min()) * 0.5
    # mid_z = (output['upper']['z'].max()+output['lower']['z'].min()) * 0.5
    # ax.set_xlim(mid_x - max_range, mid_x + max_range)
    # ax.set_ylim(mid_z - max_range, mid_z + max_range)
    # ax.set_zlim(mid_y - max_range, mid_y + max_range)
    plt.xlabel('x')
    plt.ylabel('z')

    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    
    for surface in output:
        x,y,z = output[surface].T
        ax.scatter(x, z)

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