from scipy.interpolate import interp1d
import numpy as np
import math
import pickle

from aeropy.geometry.airfoil import CST, create_x
from aeropy.geometry.wing import Rotation_euler_vector
from aeropy.CST_3D.meshing import uniform_mesh_generator
from aeropy.filehandling.vtk import generate_surface

class ControlPoints():
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

class CST_Object():
    def __init__(self, B = {'upper':[[1]]}, cp=ControlPoints(), origin = [0,0,0], 
                    axis_order = [0,1,2], mesh = (10,10)):
        self.surfaces = list(B.keys())
        self.B = B
        self.cp = cp
        self.origin = origin
        self.axis_order = axis_order
        
    def calculate_surface(self, mesh=(10,10), mesh_type='uniform'):
        if type(mesh) == dict or mesh_type=='uniform':
            self.mesh_input = mesh
            self.dimensions = list(mesh)
        else:
            self.mesh_input = {}
            for surface in self.surfaces:
                self.mesh_input[surface] = np.array(mesh)
            #TODO: how to determine dimensions for cases where I do
            # not how the mesh is (seed per x and y)
            self.dimensions = NotImplementedError
        if mesh_type is None:
            self.mesh_surface = CST_3D(self.B, self.mesh_input, 
                                       cp = self.cp, origin = self.origin, 
                                       axis_order = self.axis_order)
        else:
            self.mesh_surface = CST_3D(self.B, self.mesh_input, 
                                       cp = self.cp, origin = self.origin, 
                                       axis_order = self.axis_order, 
                                       mesh_type = mesh_type)
    
    def generate_vtk(self, filename = 'test'):
        # Format data
        data = np.reshape(self.mesh_surface, self.dimensions + [3])
        # Generate vtk
        generate_surface(data=data, filename=filename) 

        
def CST_3D(B, mesh = (10,10), cp = ControlPoints(),
           origin = [0,0,0], axis_order = [0,1,2], **optional_arguments):
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
    - optional_arguments:
        - format: controls format of output as either 'array' (size = [N,3])
                  or 'matrix' (size = [N,M,3])
        - meshtype: if equal to:
            - 'uniform' use list (len=2 or len=4) as the seeding
            - 'parameterized' uses point as is as mesh
            - 'assembly'
            - 'mixed'
    
    """
    # Process mesh
    mesh, dimensions = process_mesh(mesh, cp, origin, axis_order,
                                    optional_arguments)

    output = {}
    for surface in B:
        # Calculating non-dimensional surface
        data_parm = parameterized_transformation(mesh[surface], 
                                 B=B[surface], cp=cp, surface=surface)

        # Dimensionalize
        data_phys= dimensionalize(data_parm, cp)
        
        # Position part in assembly_transformation
        data_assembly = assembly_transformation(data_phys, axis_order, origin)
        
        # Reformat if necessary
        try:
            if optional_arguments['format'] == 'matrix':
                data_assembly.resize(dimensions[surface] + [3])
        except:
            pass
        output[surface] = data_assembly

    if len(B.keys()) == 1:
        return output[list(B.keys())[0]]
    else:
        return output

def parameterized_transformation(mesh, B, cp=ControlPoints(),
                                 surface = 'upper'):
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
        
def dimensionalize(Raw_data, cp = ControlPoints(), inverse = False):
    '''Converts from parameterized domain to physical domain'''
        
    if inverse:
        data = physical_transformation(Raw_data, cp, inverse)
        data = intermediate_transformation(data, cp, inverse)
    
    else:
        data = intermediate_transformation(Raw_data, cp, inverse)
        data = physical_transformation(data, cp, inverse)

    return data
    

def intermediate_transformation(data_nondimensional, cp = ControlPoints(), 
                                inverse = False):
    '''Converts from Parameterized domain to intermediate domain (not rotated).
       - psi, xi, eta are the coordinates of the points to be analyzed
       - eta_control points where we know values such as chord etc'''

    psi, eta, zeta = data_nondimensional.T
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

def physical_transformation(data, cp = ControlPoints(), inverse = False):
            
    if inverse:
        # CURRENTLY ONLY WORKS FOR NO TWIST ON FUSELAGE
        eta= data[:,1]/cp.half_span

        #Move to rotation and remove shear
        data[:,0] += cp.twist_origin['x']

        # Adjust if only two coordinates are provided. Using zero values
        # which is true if there is no twist
        try:
            data[:,2] += cp.twist_origin['z'] - cp.fshear(eta)
        except:
            dummy_z = np.reshape(np.zeros(data.shape[0]),(data.shape[0],1))
            data = np.append(data,dummy_z, axis=1)

        # Twisting
        R = np.transpose(Rotation_euler_vector(cp.twist_axis, cp.ftwist(eta)),(1,0,2))
        data = np.array([np.dot(R[:,:,l],data[l]) for l in range(len(eta))])

        # Returning from twist origin
        data[:,0] += cp.twist_origin['x']
        data[:,2] += cp.twist_origin['z']

    else:
        eta= data[:,1]/cp.half_span
        #displacing to center of rotation
        data[:,0] -= cp.twist_origin['x']
        data[:,2] -= cp.twist_origin['z']

        # Twisting

        R = Rotation_euler_vector(cp.twist_axis, cp.ftwist(eta))
        data = np.array([np.dot(R[:,:,l],data[l]) for l in range(len(eta))])

        #Return to position and shear
        data[:,0] -= cp.twist_origin['x']
        data[:,2] += -cp.twist_origin['z'] + cp.fshear(eta)
    return(data)

def assembly_transformation(data, axis_order, origin, inverse=False):
    '''Transform from physical domain to assembly domain'''
    if inverse:
        data -= np.array(origin)
        data[:,[0,1,2]] = data[:,axis_order]

    else:
        data[:,axis_order] = data[:,[0,1,2]]
        data += np.array(origin)
            
    return(data)
       
def inverse_problem(coordinates = {'psi':0.5, 'zeta':0.5}, B = {'upper':[[1]]}, 
                    cp=ControlPoints()):
    '''
    TODO: Inverse problem function is still NOT robust to
    function for any case. For example:
        Au = np.array([0.172802, 0.167353, 0.130747, 0.172053, 0.112797, 0.168891])
        Al = np.array([0.163339, 0.175407, 0.134176, 0.152834, 0.133240, 0.161677])
        # Wing properties (control points)
        cp = ControlPoints()
        cp.set(chord = [1.,1.],
               twist = [0, 0])'''
               
    from scipy.optimize import minimize
    
    def residual(psi, eta, zeta):
        mesh = {surface: np.array([[psi[0],eta[0]],])}
        output = CST_3D(B, mesh = mesh, cp = ControlPoints(),
                        mesh_type = 'parameterized')  
        return abs(zeta - output[0][2])
    
    zeta = coordinates['zeta']
    surface = list(B.keys())[0]
    if 'psi' in coordinates:
        psi = coordinates['psi']
        f = lambda x: residual([psi],x,zeta)
        solution = minimize(f,0.5, bounds = [[0,1]])
        eta = solution.x[0]
    elif 'eta' in coordinates:
        eta = coordinates['eta']
        f = lambda x: residual(x,[eta],zeta)
        solution = minimize(f,0.5, bounds = [[0,1]])
        psi = solution.x[0]

    # Check if found actual solution
    tol = 1e-4

    if abs(solution.fun) < tol:
        return psi, eta, zeta
    else:
        raise Exception('There is no point on surface with required z coordinate')


def process_mesh(mesh, cp = ControlPoints(), origin = [0,0,0],
                 axis_order = [], options={}):
    ''' 'Universal' tool processing for following cases:
        - mesh_options not defined: mesh=[Npsi_upper, Neta_upper, Npsi_lower, Neta_lower]
                                    or mesh=[Npsi_upper, Neta_upper]=[Npsi_lower, Neta_lower]
                                    '''

    if options['mesh_type'] == 'uniform':
        mesh, dimensions = uniform_mesh_generator(mesh)
    else:
        if options['mesh_type'] == 'parameterized':
            dimensions = NotImplementedError
            
        elif options['mesh_type'] == 'mixed':
            # non-dimensionalize y
            for surface in mesh:
                mesh[surface][:,1] /= cp.half_span
            dimensions = NotImplementedError
            
        elif options['mesh_type'] == 'physical':
            for surface in mesh:
                mesh[surface] = dimensionalize(mesh[surface], cp = cp, 
                                               inverse = True)[:,[0,1]]
            dimensions = NotImplementedError
        
        elif options['mesh_type'] == 'assembly':
            for surface in mesh:
                mesh[surface] = assembly_transformation(mesh[surface], axis_order, origin, inverse=True)
                mesh[surface] = dimensionalize(mesh[surface], cp = cp, 
                                               inverse = True)[:,[0,1]]
            dimensions = NotImplementedError

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
    cp = ControlPoints()
    cp.set(chord = [1.,.3],
           twist = [0, .3],
           shear = [0, .1],
           sweep = [0, .3],
           N1 = [.5, 1.],
           N2 = [1., 1.])
 
    B = {'upper':[Au,Au], 'lower':[Al,Al]}
    
    # Using CST_3D function
    #output = CST_3D(B, cp = cp, format='matrix')
    #generate_vtk(output, filename='test')
    
    # Using CST_object
    wing_upper = CST_Object(B = {'upper':[Au,Au]}, cp = cp)
    wing_upper.calculate_surface((20,20))
    wing_upper.generate_vtk('upper')
    
    wing_lower = CST_Object(B = {'lower':[Al,Al]}, cp = cp)
    wing_lower.calculate_surface((20,20))
    wing_lower.generate_vtk('lower')
    
    # Just for plotting
    output = {'upper':wing_upper.mesh_surface,
              'lower':wing_lower.mesh_surface}
    
    # Inverse problem
    # eta = 0.5
    # psi_list = []
    # eta_list = []
    # zeta_list = []
    # for zeta in np.linspace(0,0.03, 10):
        # [psi, eta, zeta] = inverse_problem(coordinates = {'eta':eta, 'zeta':zeta}, 
                                           # B = {'upper':[Au,Au]}, 
                                           # cp=cp)

        # psi_list.append(psi)
        # eta_list.append(eta)
        # zeta_list.append(zeta)
    
    # Plotting
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # # ax.scatter(psi_list,  zeta_list, eta_list)
    # for surface in output:
        # x,y,z = output[surface].T
        # ax.scatter(x, z, y, cmap=plt.get_cmap('jet'),
                        # linewidth=0, antialiased=False)
    # ax.set_zlim(0, 4)
    # plt.xlabel('x')
    # plt.ylabel('z')

    # fig2 = plt.figure()
    # ax = fig2.add_subplot(111)
    # for surface in output:
        # x,y,z = output[surface].T
        # ax.scatter(x, z)
    # plt.xlabel('x')
    # plt.ylabel('z')
    # plt.show()
