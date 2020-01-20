import math
import numpy as np
from scipy import integrate, optimize


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
