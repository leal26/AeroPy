import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize

class polynomial():
    """class for a polynomial function
    """

    def __init__(self, D=[0, -math.sqrt(3)/3, math.sqrt(3)/3, 0], chord = None,
                 color = 'k'):
        self.D = D
        self.chord = chord
        self.color = color

    def x1(self, x1, diff=None):
        if diff is None:
            return x1
        elif diff == 'x1':
            return np.ones(len(x1))
        elif diff == 'x11':
            return np.zeros(len(x1))
        elif diff == 'theta3' or diff == 'theta33':
            return np.zeros(len(x1))
        elif diff == 'theta1':
            dr = self.r(x1, diff='x1')
            return 1/np.sqrt(np.einsum('ij,ij->i',dr, dr))
        elif diff == 'theta11':
            return -self.x1(x1, 'theta1')**4*self.x3(x1, 'x1')*self.x3(x1, 'x11')

    def x2(self, x1, diff=None):
        return np.zeros(len(x1))

    def x3(self, x1, diff=None, D=None):
        """ z2 (checked)"""
        if D is None:
            D = self.D

        if len(D) != 5:
            D_temp = np.zeros(5)
            D_temp[:len(D)] = D
            D = D_temp

        if diff is None:
            return(D[4]*x1**4 + D[3]*x1**3 + D[2]*x1**2 + D[1]*x1 + D[0])
        elif diff == 'x1':
            return(4*D[4]*x1**3 + 3*D[3]*x1**2 + 2*D[2]*x1 + D[1])
        elif diff == 'x11':
            return(12*D[4]*x1**2 + 6*D[3]*x1 + 2*D[2])
        elif diff == 'x111':
            return(24*D[4]*x1 + 6*D[3])
        elif diff == 'theta3':
            return(np.zeros(len(x1)))

    def r(self, x1 = None, diff=None):
        if x1 is None:
            x1 = self.x1_grid
        else:
            if type(x1) == float:
                x1 = np.array([x1])
        if diff == 'theta1':
            output = np.array([self.x1(x1, 'x1'),
                               self.x2(x1, 'x1'),
                               self.x3(x1, 'x1')]).T
            output = np.einsum('ij,i->ij',output, self.x1(x1, 'theta1'))
        else:
            output = np.array([self.x1(x1, diff),
                               self.x2(x1, diff),
                               self.x3(x1, diff)]).T
            self.position = output
        return (output)

    def basis(self, x1 = None, diff=None):
        if x1 is None:
            x1 = self.x1_grid

        if diff is None:
            self.a = np.zeros([3,  len(x1), 3])
            self.a[0,:,:] = np.array([self.x1(x1, 'x1')*self.x1(x1, 'theta1'),
                                     [0]*len(x1),
                                      self.x3(x1, 'x1')*self.x1(x1, 'theta1')]).T
            self.a[1,:,:]  = np.array([[0,1,0],]*len(x1))
            self.a[2,:,:]  = np.cross(self.a[0,:,:], self.a[1,:,:])

        elif diff == 'theta':
            # Most components are null
            self.da = np.zeros([3,3,len(x1),3])
            #  a1 diff theta1
            self.da[0,0,:,:] = np.einsum('ij,i->ij', self.r(x1, 'x11'), self.x1(x1, 'theta1')**2) + \
                               np.einsum('ij,i->ij', self.r(x1, 'theta1'), self.x1(x1, 'theta11'))
            #  a3 diff theta1
            self.da[0,2,:,:] = np.einsum('ij,i->ij', self.r(x1, 'x11'), self.x1(x1, 'theta1')**2) + \
                               np.einsum('ij,i->ij', self.r(x1, 'theta3'), self.x1(x1, 'theta33'))

    def christoffel(self, i, j, k, order=1):
        if order == 1:
            gik_j = self.dA[i,k,j]
            gjk_i = self.dA[j,k,i]
            gij_k = self.dA[i,j,k]
            return .5*(gik_j + gjk_i - gij_k)
        elif order == 2:
            raise NotImplementedError

    def metric_tensor(self, diff = None):

        if diff is None:
            self.A = np.zeros([3,3,len(self.x1_grid)])
            for i in range(3):
                for j in range(3):
                    self.A[i,j] = np.einsum('ij,ij->i',self.a[i,:], self.a[j,:])

        elif diff == 'theta':
            self.dA = np.zeros([3,3,3,len(self.x1_grid)])
            for i in range(3):
                for j in range(3):
                        for k  in range(3):
                            self.dA[i,j,k] = np.einsum('ij,ij->i',self.da[i,k,:],
                                                        self.a[j,:]) + \
                                             np.einsum('ij,ij->i',self.a[i,:],
                                                        self.da[j,k,:])

    def curvature_tensor(self):
        self.B = np.zeros([2,2,len(self.x1_grid)])
        for alpha in range(2):
            for beta in range(2):
                self.B[alpha, beta] = self.christoffel(alpha, beta, 2)

    def arclength(self, chord = None):
        def integrand(x1):
            dr = self.r(x1, 'x1')
            return np.sqrt(np.inner(dr, dr))
        if chord is None:
            chord = self.chord
        return integrate.quad(integrand, 0, chord)

    def calculate_x1(self, length_target):
        def f(c_c):
            length_current, err = self.arclength(c_c)
            return abs(target - length_current)

        x1 = []
        for target in length_target:
            x1.append(optimize.minimize(f, target).x[0])
        self.x1_grid = np.array(x1)

    def plot(self, basis=False, label=None, linestyle = '-', color = None):
        r = self.r(self.x1_grid)

        if color is None:
            color = self.color

        if label is None:
            plt.plot(r[:,0], r[:,2], color, linestyle = linestyle, lw = 4)
        else:
            plt.plot(r[:,0], r[:,2], color, linestyle = linestyle, lw = 4,
                     label=label)
        if basis:
            plt.quiver(r[:,0], r[:,2],
                       self.a[0,:,0], self.a[0,:,2],
                       angles='xy', color = self.color, scale_units='xy')
            plt.quiver(r[:,0], r[:,2],
                       self.a[2,:,0], self.a[2,:,2],
                       angles='xy', color = self.color, scale_units='xy')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')

class poly():
    """class for a polynomial function
    """

    def __init__(self, a=[0, -math.sqrt(3)/3, math.sqrt(3)/3, 0]):
        self.a = a

    def z2(self, z1, diff=None, a=None):
        """ z2 (checked)"""
        if a is None:
            a = self.a
        if diff is None:
            return(a[3]*z1**3 + a[2]*z1**2 + a[1]*z1 + a[0])
        elif diff == 'z1':
            return(3*a[3]*z1**2 + 2*a[2]*z1 + a[1])
        elif diff == 'z11':
            return(6*a[3]*z1 + 2*a[2])
        elif diff == 'z111':
            return(6*a[3])
        elif diff == 'x1':
            return(self.z2(z1, 'z1')*self.z1(z1, 'x1'))
        elif diff == 'x11':
            return(self.z2(z1, 'z11')*(self.z1(z1, 'x1'))**2 +
                   self.z2(z1, 'z1')*self.z1(z1, 'x11'))

    def x1(self, z1, diff=None, a=None):
        """ dx1/ dz1 (checked)"""
        if diff is None:
            output = []
            try:
                for z_final in z1:
                    output_i, err = integrate.quad(lambda x: self.x1(x, 'z1'),
                                                   0, z_final)
                    output.append(output_i)
                output = np.array(output)
            except(TypeError):
                output, err = integrate.quad(lambda x: self.x1(x, 'z1'), 0, z1)
            return(output)
        elif diff == 'z1':
            return(np.sqrt(1+(self.z2(z1, 'z1', a=a))**2))
        elif diff == 'z11':
            return(self.z1(z1, 'x1')*self.z2(z1, 'z1')*self.z2(z1, 'z11'))

    def z1(self, input, diff=None, a=None):
        """ dx1 / dz1 (all checked). For calculating z1 from x1, there is not
           a numerical solution, but I can minimize the residual."""
        if diff is None:
            output = []
            for x_final in input:
                def _residual(x):
                    return abs(x_final - self.x1(x))
                output_i = optimize.newton(_residual, x_final)
                output.append(output_i)
            output = np.array(output)
            return(output)
        if diff == 'x1':
            if a is None:
                return(1.0/self.x1(input, 'z1'))
            else:
                return(- self.z1(input, 'z1')**3*self.z2(input, 'z1') *
                       self.z2(input, 'z1', a=a))
        elif diff == 'x11':
            return(-self.z1(input, 'x1')**4*self.z2(input, 'z1') *
                   self.z2(input, 'z11'))
        elif diff == 'x111':
            return(-4*self.z1(input, 'x1')**3*self.z1(input, 'x11') *
                   self.z2(input, 'z1') * self.z2(input, 'z11') -
                   self.z1(input, 'x1')**5*(self.z2(input, 'z11')**2 -
                                            self.z2(input, 'z1') *
                                            self.z2(input, 'z111')))

    def tangent(self, z1, diff=None):
        """Tangent vector r (checked)"""
        if diff is None:
            try:
                output = self.z1(z1, 'x1')*np.array([np.ones(len(z1)),
                                                     self.z2(z1, 'z1')])
            except(ValueError):
                output = self.z1(z1, 'x1')*np.array([1, self.z2(z1, 'z1')])
            if len(output) > 2:
                BRAKE
        elif diff == 'x1':
            try:
                output = self.z1(z1, 'x1')**2*np.array([np.zeros(len(z1)),
                                                        self.z2(z1, 'z11')])
                output += self.z1(z1, 'x11')*np.array([np.ones(len(z1)),
                                                       self.z2(z1, 'z1')])
            except(ValueError):
                output = self.z1(z1, 'x1')**2*np.array([0, self.z2(z1, 'z11')])
                output += self.z1(z1, 'x11')*np.array([1, self.z2(z1, 'z1')])

        return(output)

    def normal(self, z1, diff=None):
        """Normal vector (checked)"""
        if diff is None:
            try:
                output = self.z1(z1, 'x1')*np.array([- self.z2(z1, 'z1'),
                                                     np.ones(len(z1))])
            except(TypeError):
                output = self.z1(z1, 'x1')*np.array([- self.z2(z1, 'z1'), 1])
        elif diff == 'z1':
            try:
                output = self.z1(z1, 'x1')*np.array([- self.z2(z1, 'z1'),
                                                     np.zeros(len(z1))])
            except(TypeError):
                output = self.z1(z1, 'x1')*np.array([- self.z2(z1, 'z11'), 0])
        elif diff == 'x1':
            output = self.z1(z1, 'x11')*np.array([- self.z2(z1, 'z1'),
                                                  np.ones(len(z1))])
            output += self.z1(z1, 'x1')*self.normal(z1, 'z1')
        elif diff == 'x11':
            output = self.z1(z1, 'x111')*np.array([- self.z2(z1, 'z1'),
                                                   np.ones(len(z1))])
            output += self.z1(z1, 'x11')*np.array([- self.z2(z1, 'z11'),
                                                   np.zeros(len(z1))])
            output += self.z1(z1, 'x11')*self.normal(z1, 'z1')
            output += self.z1(z1, 'x1')**2*self.normal(z1, 'z1')
        elif type(diff) == list:
            output = self.z1(z1, 'x1')*np.array([- self.z2(z1, 'z1', a=diff),
                                                 np.zeros(len(z1))])
            output += self.z1(z1, 'x1', a=diff)*np.array([- self.z2(z1, 'z1'),
                                                          np.ones(len(z1))])
        elif diff == 'x2':
            output = np.array([[0], [0]])
        return(output)

    def neutral_line(self, z1, a=None):
        """ Position along neutral line"""
        if a is not None:
            return(np.array([z1, self.z2(z1, a=a)]))
        else:
            return(np.array([z1, self.z2(z1)]))

    def r(self, input, x2=0, diff=None, input_type='x1'):
        """ Position anywhere along shell considering shell thickness """

        z1 = self._process_input(input, input_type)
        if diff is None:
            output = self.neutral_line(z1) + x2*self.g(2, z1)
        elif diff == 'x1':
            output = self.g(1, z1, x2)
        elif diff == 'x2':
            output = self.g(2, z1, x2)
        elif type(diff) == list:
            output = self.neutral_line(z1, a=diff) + \
                x2*self.normal(z1, diff=diff)
        return(output)

    def g(self, i, input, x2=0, diff=None, input_type='z1'):
        """Tangent vector r (checked)"""
        z1 = self._process_input(input, input_type)

        if i == 1:
            if diff is None:
                g_i = self.tangent(z1) + x2*self.normal(z1, diff='x1')
            elif diff == 'x1':
                g_i = self.tangent(z1, 'x1') + x2*self.normal(z1, diff='x2')
            elif diff == 'x2':
                g_i = self.normal(z1, diff='x1')
        elif i == 2:
            if diff is None:
                g_i = self.normal(z1)
            else:
                g_i = np.array([[0], [0]])
        return(g_i)

    def gij(self, i, j, z1, x2=0, diff=None, covariant=True, orthogonal=True):
        def dot(a, b):
            a_1, a_2 = a
            b_1, b_2 = b
            return(a_1*b_1 + a_2*b_2)

        if diff is None:
            gi = self.g(i, z1, x2)
            gj = self.g(j, z1, x2)
            output = []
            for n in range(len(z1)):
                output.append(gi[0][n]*gj[0][n] +
                              gi[1][n]*gj[1][n])

            output = np.array(output)

        elif diff == 'x1':
            # Calculating basic vectors
            t = self.tangent(z1)
            dt = self.tangent(z1, 'x1')
            n = self.normal(z1)
            dn = self.normal(z1, 'x1')
            ddn = self.normal(z1, 'x11')
            # calculate g
            if i == 1 and j == 1:
                output = 2*dot(t, dt) + 2*x2*(dot(dt, dn) + dot(t, ddn)) + \
                    2*x2**2*dot(dn, ddn)
            if (i == 1 and j == 2) or (i == 2 and j == 1):
                output = x2*(dot(dn, dn) + dot(n, ddn))
            if i == 2 and j == 2:
                output = 2*dot(n, dn)
        elif diff == 'x2':
            # Calculating basic vectors
            t = self.tangent(z1)
            n = self.normal(z1)
            dn = self.normal(z1, 'x1')
            # calculate g
            if i == 1 and j == 1:
                output = 2*dot(t, dn) + 2*x2**2*dot(dn, dn)
            if (i == 1 and j == 2) or (i == 2 and j == 1):
                output = dot(dn, n)
            if i == 2 and j == 2:
                output = np.zeros(len(z1))
        if not covariant and orthogonal:
            if i == j:
                output = 1/output
            else:
                output = 0.0
        return(output)

    def christoffel(self, i, k, l, z1, x2, order='second'):
        if order == 'second':
            output = np.zeros(len(z1))
            for m in range(1, 3):
                gim = self.gij(i, m, z1, x2, covariant=False)
                gmk_l = self.gij(m, k, z1, x2, diff='x%i' % (l))
                gml_k = self.gij(m, l, z1, x2, diff='x%i' % (k))
                gkl_m = self.gij(k, l, z1, x2, diff='x%i' % (m))
                output += .5*gim*(gmk_l + gml_k + gkl_m)
        return(output)

    def _process_input(self, input, input_type):
        if input_type == 'z1':
            z1 = input
        elif input_type == 'x1':
            z1 = self.z1(input)
        return(z1)
