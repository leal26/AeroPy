import aeropy
import numpy as np
from scipy import optimize

class shell():
    def __init__(self, geometry_parent, geometry_child,  properties,
                 bc, chord=1, ndim=2):
        """
        COnstant length is assumed for the structure

        A: metric tensor
        dA: metric tensor covariant derivative as a function of theta
        a: curvilinear basis vectors
        chord: length of the beam"""
        # Defining geometries
        self.g_p = geometry_parent
        self.g_c = geometry_child

        # shell thickness
        self.h = properties.dimensions[1]
        self.width = properties.dimensions[0]
        self.ndim = 2
        self.bc = bc
        self.properties = properties

    def kinematics_p(self, x1):
        self.g_p.basis(x1)

    def kinematics_c(self, x1):
        # define new chord
        chord = aeropy.CST_2D.calculate_c_baseline(c_L, Au_C, Au_L, deltaz)
        x1 = chord*x1
        # calculate basis vectors
        self.g_c.basis(x1)

    def calculate_chord(self, length_target = None):
        def f(c_c):
            length_current, err = self.g_c.arclength(c_c)
            return abs(length_target - length_current)
        if length_target is None:
            length_target, err = self.g_p.arclength()
        self.g_c.chord = optimize.minimize(f, self.g_p.chord).x[0]
        # In case the calculated chord is really close to the original
        if abs(self.g_p.chord - self.g_c.chord) < 1e-7:
            self.g_c.chord = self.g_p.chord

    def calculate_strains(self):
        self.gamma = 0.5*(self.g_c.A[:2, :2, :] - self.g_p.A[:2, :2, :])

    def calculate_change_curvature(self):
        self.rho = -(self.g_c.B - self.g_p.B)

    def CauchyGreen(self):
        """From the definition of Hookian Thin homogeneous isentropic shell
        (Eq. 9.98a) from Wempner's book:
            - the definition uses contravariant basis vectors, but this whole
              class uses covariant basis vectors. because of that all values
              are inverted (coordinate system assumed orthogonal)"""
        self.C = np.zeros([2,2,2,2,len(self.g_c.x1_grid)])
        c0 = self.properties.young/2/(1+self.properties.poisson)
        for alpha in range(2):
            for beta in range(2):
                for gamma in range(2):
                    for eta in range(2):
                        a1 = self.g_c.A[alpha, gamma]*self.g_c.A[beta, eta]
                        a2 = self.g_c.A[alpha, eta]*self.g_c.A[beta, gamma]
                        a3 = self.g_c.A[alpha, beta]*self.g_c.A[gamma, eta]
                        c3 = (2*self.properties.poisson)/(1-self.properties.poisson)
                        self.C[alpha, beta, gamma, eta, :] = c0*(a1 + a2 + c3*a3)

    def free_energy(self):
        self.phi_M = (self.h/2)*np.einsum('ijklm,ijm,klm->m',self.C,self.gamma,self.gamma)
        self.phi_B = (self.h**3/24)*np.einsum('ijklm,ijm,klm->m',self.C,self.rho,self.rho)

        self.phi =  self.phi_B + self.phi_M


    def strain_energy(self):
        # print('M', self.phi_M)
        # print('B', self.phi_B)
        self.U = self.width*np.trapz(self.phi, self.g_p.x1_grid)

    def work(self):
        energy = 0
        for i in range(self.bc.concentrated_n):
            u = self.g_c.position[-1,:] - self.g_p.position[-1,:]
            print('c',self.g_c.position[-1,:] )
            print('p',self.g_p.position[-1,:] )
            print('u', u)
            print(self.bc.concentrated_load)
            for j in range(3):
                energy += self.bc.concentrated_load[i][j] * u[j]
        self.W = energy

    def residual(self):
        self.R = abs(self.U - self.W)

class design_exploration():
    def __init__(self, component):
        pass

    def sweep_geometries(self, geom_variables, input_function, reorder=None,
                         loading_condition = 'plane_stress'):
        energy_list = []
        residual_list = []
        n = len(geom_variables)
        for i in range(n):
            print(i)
            input = geom_variables[i]
            residual_list.append(self.residual(input, input_type = 'Geometry',
                                               input_function = input_function,
                                               loading_condition = loading_condition))
            energy_list.append(self.strain_energy())
        if reorder is not None:
            residual_list = np.resize(residual_list, reorder)
            energy_list = np.resize(energy_list, reorder)
        return(energy_list, residual_list)
