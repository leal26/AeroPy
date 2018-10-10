import math
import numpy as np
from scipy import integrate, optimize


class poly():
    """class for a polynomial function
    """

    def __init__(self, a=[0, -math.sqrt(3)/3, math.sqrt(3)/3, 0], alpha=1,
                 config='parent'):
        self.a = a
        self.alpha = alpha
        self.config = config

    def z2(self, z1, diff=None, a=None):
        """ z2 (checked)"""
        if a is None:
            a = self.a
        if diff is None:
            return(a[0]*z1**3 + a[1]*z1**2 + a[2]*z1 + a[3])
        elif diff == 'z1':
            return(3*a[0]*z1**2 + 2*a[1]*z1 + a[2])
        elif diff == 'z11':
            return(6*a[0]*z1 + 2*a[1])
        elif diff == 'z111':
            return(6*a[0])
        elif diff == 'x1':
            return(self.z2(z1, 'z1')*self.z1(z1, 'x1'))
        elif diff == 'x11':
            return(self.z2(z1, 'z11')*(self.z1(z1, 'x1'))**2 +
                   self.z2(z1, 'z1')*self.z1(z1, 'x11'))

    def inflection_points(self):
        return([-self.a[1]/(3*self.a[0])])

    def x1(self, z1, diff=None):
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
            return(np.sqrt(1+(self.z2(z1, 'z1'))**2))
        elif diff == 'z11':
            return(self.z1(z1, 'x1')*self.z2(z1, 'z1')*self.z2(z1, 'z11'))

    def z1(self, input, diff=None):
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
            return(1.0/self.x1(input, 'z1'))
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
        return(output)

    def neutral_line(self, z1):
        """ Position along neutral line"""
        return(np.array([z1, self.z2(z1)]))

    def r(self, input, x2=0, normalize=False, diff=None, input_type='x1'):
        """ Position anywhere along shell considering shell thickness """

        if input_type == 'z1':
            z1 = input
        elif input_type == 'x1':
            if self.config == 'parent':
                z1 = self.z1(input)
            elif self.config == 'child':
                z1 = self.z1(self.alpha*input)

        if diff is None:
            output = self.neutral_line(z1) + x2*self.g(2, z1)
        elif diff == 'x1':
            output = self.alpha*self.g(1, z1, x2)
        elif diff == 'x2':
            output = self.g(2, z1, x2)
        elif type(diff) == list:
            a = self.a.copy()
            for i in range(len(diff)):
                if diff[i] == 0:
                    a[i] = 0.0
            output = np.array([z1, self.z2(z1, a=a)])
        if normalize:
            output = output.T
            for i in range(len(output)):
                output[i] /= np.linalg.norm(output[i])
            output = output.T
        return(output)

    def g(self, i, z1, x2=0, diff=None):
        """Tangent vector r (checked)"""
        if i == 1:
            g_i = self.tangent(z1) + x2*self.normal(z1, diff='x1')
        elif i == 2:
            g_i = self.normal(z1)
        if diff is None:
            return(g_i)
        elif diff == 'x1':
            return

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


class frame():
    def __init__(self, curve=poly(), frame='Frenet-Serret',
                 z1=np.linspace(0, 1, 11)):
        self.curve = poly()
        self.frame = frame
        self.z1 = z1
        self.z2 = self.curve.z2(z1)

    def z1_default(self, z1):
        if z1 is None:
            return(self.z1)
        else:
            return(z1)

    def bishop(self, z1=None):
        z1 = self.z1_default(z1)
        T = self.T(z1)
        alpha = self.alpha(z1)

        M1 = self.M1(z1=z1, alpha=alpha)
        M2 = self.M2(z1=z1, alpha=alpha)

        return(T, M1, M2)

    def frenet_serret(self, z1=None):
        z1 = self.z1_default(z1)
        T = self.T(z1)
        N = self.N(z1=z1)
        B = self.B(z1=z1)

        return(T, N, B)

    def T(self, z1=None):
        z1 = self.z1_default(z1)
        return(self.curve.rx1(z1))

    def N(self, z1=None):
        z1 = self.z1_default(z1)
        return(self.curve.rx11(z1))

    def M1(self, z1=None, alpha=None):
        z1 = self.z1_default(z1)
        if alpha is None:
            alpha = self.alpha(z1)
        if self.frame == 'Bishop':
            new = np.cos(alpha)*self.N(z1) - np.sin(alpha)*self.B(z1)
            return(np.cos(alpha)*self.N(z1) -
                   np.sin(alpha)*self.B(z1))

    def M2(self, z1=None, alpha=None):
        z1 = self.z1_default(z1)
        if alpha is None:
            alpha = self.alpha(z1)
        if self.frame == 'Bishop':
            return(np.sin(alpha)*self.N(z1) +
                   np.cos(alpha)*self.B(z1))

    def B(self, z1=None):
        z1 = self.z1_default(z1)
        try:
            return(np.ones(len(z1)))
        except(ValueError):
            return(1)

    def curvature(self, z1=None):
        z1 = self.z1_default(z1)
        return(self.curve.rx11(z1))

    def torsion(self, z1=None):
        z1 = self.z1_default(z1)
        try:
            return(np.zeros(len(z1)))
        except(ValueError):
            return(0)

    def alpha(self, z1=None, alpha0=0.0):
        z1 = self.z1_default(z1)
        alpha_list = []
        inflections = self.curve.inflection_points()
        j = 0
        alpha_j = alpha0
        for i in range(len(z1)):
            if j != len(inflections):
                if z1[i] >= inflections[j]:
                    alpha_j += math.pi
                    j += 1
            alpha_list.append(alpha_j)
        return(np.array(alpha_list))


def B(x, k, i, t):
    if k == 0:
        return 1.0 if t[i] <= x < t[i+1] else 0.0
    if t[i+k] == t[i]:
        c1 = 0.0
    else:
        c1 = (x - t[i])/(t[i+k] - t[i]) * B(x, k-1, i, t)
    if t[i+k+1] == t[i+1]:
        c2 = 0.0
    else:
        c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * B(x, k-1, i+1, t)
    return c1 + c2


def bspline(x, t, c, k):
    n = len(t) - k - 1
    assert (n >= k+1) and (len(c) >= n)
    return sum(c[i] * B(x, k, i, t) for i in range(n))


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # Bezier Curve
    k = 2
    t = [0, 1, 2, 3, 4, 5, 6]
    c = [-1, 2, 0, -1]
    x = np.linspace(1.5, 4.5, 50)
    y_bezier = []
    for x_i in x:
        y_bezier.append(bspline(x_i, t, c, k))

    # Hicks-Henne
    plt.plot(x, y_bezier, label='Bezier')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.legend()
    plt.show()
