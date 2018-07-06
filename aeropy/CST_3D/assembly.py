from aeropy.CST_3D.module import *
from aeropy.CST_3D.meshing import *
from aeropy.filehandling.vtk import generate_points
import numpy as np

class Aircraft():
    def __init__(self, fuselage = None, 
                 wing_upper = None, 
                 wing_lower = None, 
                 tail = None):
        self.fuselage = fuselage
        self.wing_upper = wing_upper
        self.wing_lower = wing_lower
        self.tail = tail
    
        self.intersections = None
        
    def find_intersections(self, eta0, maxit=100):
        self.intersections = {}
        self.intersections['upper'] = intersection_curve(self.wing_upper, 
                                                         self.fuselage, eta0)
        self.intersections['lower'] = intersection_curve(self.wing_lower, 
                                                         self.fuselage, eta0)
        
    def fitting(self, ):
        return 0
        
    def meshing(self, ):
        return 0
        
    def run_panair(self, ):
        return 0
        
    def run_sboom(self, ):
        return 0
        
    def run_pyldb(self, ):
        return 0
        
    def generate_vtk(self, filename = 'test'):
        wing_upper.generate_vtk('assembly_upper')
        wing_lower.generate_vtk('assembly_lower')
        fuselage.generate_vtk('assembly_fuselage')
        
        if self.intersections is not None:
            generate_points(self.intersections['upper'], 'upper_intersection')
            generate_points(self.intersections['lower'], 'lower_intersection')
    def generate_stl(self, ):
        return 0 
        
def intersection_curve(wing, fuselage, eta0, debugging = False):
    '''Find coordinates along x (including all geometry of wing)
       that include fuselage.'''

    def _intersection_point(psi, eta0, wing, fuselage, 
                            debugging = False):
        '''Works if:
            - fuselage has no twist
            - wing chord smaller than fuselage length
            - wing is not rotated and fuselage is rotated 90 in x and z'''
        tol=1e-3
        error = 9999
        maxit = 100
        counter = 0

        while error > tol:
            if counter == maxit:
                break
            if counter == 0:
                wing.calculate_surface([[psi, eta0],], 'parameterized')
            else:
                wing.calculate_surface([[psi, fuselage.mesh_surface[0][1]],], 'mixed')

            fuselage.calculate_surface(wing.mesh_surface, 'assembly')

            error = abs(wing.mesh_surface[0][1] - fuselage.mesh_surface[0][1]) + \
                    abs(wing.mesh_surface[0][0] - fuselage.mesh_surface[0][0])

            counter += 1
        
        if debugging:
            return({'wing':wing.mesh_surface[0], 'fuselage':fuselage.mesh_surface[0]})
        else:
            return(wing.mesh_surface[0]) 
        
    # Define x
    psi = np.linspace(0,1,10)

    if debugging:
        solutions = {'wing': np.zeros((len(psi),3)),
                     'fuselage':np.zeros((len(psi),3))}
    else:
        solutions = np.zeros((len(psi),3))
    for i in range(len(psi)):
        solution = _intersection_point(psi[i], eta0, wing, fuselage,debugging)
        if debugging:
            for key in solution:
                solutions[key][i] = solution[key]
        else:
            solutions[i] = solution
    return solutions
    
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    # Shape coefficient for wing cross section
    Au = np.array([0.172802, 0.167353, 0.130747, 0.172053, 0.112797, 0.168891])
    Al = np.array([0.163339, 0.175407, 0.134176, 0.152834, 0.133240, 0.161677])
    B_w = {'upper':[Au,Au], 'lower':[Al,Al]}
    cp_w = ControlPoints()
    cp_w.set(chord = [1,.1])
    
    # Shape coefficient for fuselage
    B_f = {'upper':[[3.]]}
    half_span = 5.
    chord_control_f = .4
    cp_f = ControlPoints()
    cp_f.set(eta = [0, .1, 1.],
             N1 = [.5, .5, .5],
             N2 = [1., 1., 1.],
             chord = [0, chord_control_f, 0],
             sweep = [0.0,0.,chord_control_f/2],
             twist = [0., 0., 0.],
             shear = [0., 0., 0.],
             half_span = half_span)
          
    # Defining CST parts
    wing_upper = CST_Object(B = {'upper':[Au,Au]}, cp = cp_w)
    wing_lower = CST_Object(B = {'lower':[Au,Au]}, cp = cp_w)
    fuselage  =  CST_Object(B_f, axis_order=[2,0,1], cp=cp_f, 
                      origin = [-half_span*.3, 0,-chord_control_f/2.])
    
    # Assembly
    JAXA = Aircraft(fuselage   = fuselage,
                    wing_upper = wing_upper,
                    wing_lower = wing_lower)

    # Calculating intersections
    JAXA.find_intersections(eta0 = chord_control_f/2.)
    
    # Generating data for generic plot
    wing_upper.calculate_surface((10,10))
    wing_lower.calculate_surface((10,10))
    fuselage.calculate_surface((40,40))
    output_w = {'upper': wing_upper.mesh_surface,
                'lower': wing_lower.mesh_surface}
    output_f = fuselage.mesh_surface 

    # Generating vtk files
    JAXA.generate_vtk()
    
    # Plotting
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    for surface in output_w:
        x,y,z = output_w[surface].T
        ax.scatter(x, y, z, color='b', linewidth=0, antialiased=False)
        
        x,y,z = JAXA.intersections[surface].T
        ax.plot(x, y, z, 'g', lw=4)
        ax.scatter(x, y, z, c='g')
        
    x,y,z = output_f.T
    ax.scatter(x, y, z, color='k', linewidth=0, antialiased=False)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()