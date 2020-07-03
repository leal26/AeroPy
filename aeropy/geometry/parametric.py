import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize

from aeropy.geometry.airfoil import CST
from aeropy.CST_2D import dxi_u, ddxi_u


class CoordinateSystem(object):
    """class for a polynomial function
    """

    def __init__(self, D, **kwargs):
        # Other inputs vary according to shape class
        for key, value in kwargs.items():
            setattr(self, key, value)
        # D is always the shape parameters
        self.D = D

    @property
    def D(self):
        return self._D

    @D.setter
    def D(self, values):
        if self.n == 1 and (isinstance(values, float) or isinstance(values, np.float64)):
            self._D = values
        else:
            print(values, self.n, len(values))
            if len(values) != self.n:
                self._D = np.zeros(self.n)
                self._D[:len(values)] = values
            else:
                self._D = values

    def x1(self, x1, diff=None):
        if diff is None:
            return x1
        elif diff == 'x1':
            try:
                return np.ones(len(x1))
            except:
                return 1
        elif diff == 'x11':
            try:
                return np.zeros(len(x1))
            except:
                return 0
        elif diff == 'theta3':
            return(self.a[2, :, 0])
        elif diff == 'theta1':
            output = 1/np.sqrt(1 + self.x3(x1, 'x1')**2)
            output[np.isnan(output)] = 0
            return output
            # return np.ones(len(x1))
        elif diff == 'theta11':
            # return -self.x1(x1, 'theta1')**4*self.x3(x1, 'x1')*self.x3(x1, 'x11')
            # dr = self.r(x1, diff='x1')
            # ddr = self.r(x1, diff='x11')
            # a1 = 1/np.sqrt(np.einsum('ij,ij->i',dr, dr))**3
            # a2 = np.einsum('ij,ij->i',dr, ddr)
            # return np.multiply(a1, a2)*self.x1(x1, 'theta1')
            # return np.zeros(len(x1))
            return -self.x1(x1, 'theta1')**4*self.x3(x1, 'x1')*self.x3(x1, 'x11')
        elif diff == 'theta31' or diff == 'theta13':
            return(-self.x3(x1, 'theta11'))
            # return(np.zeros(len(x1)))
        else:
            return(np.zeros(len(x1)))

    def x2(self, x1, x2_value=0, diff=None):
        if diff is None:
            return np.zeros(len(x1))
        if diff is 'x2' or diff is 'theta2':
            return np.ones(len(x1))

    @classmethod
    def polynomial(cls, D, chord=1, color='b', n=6):
        c = cls(D, chord=chord, color=color, n=n)
        c.x3 = c._x3_poly
        return c

    @classmethod
    def CST(cls, D, chord, color, N1=.5, N2=1.0, deltaz=0):
        c = cls(D, chord=chord, color=color, n=len(D), N1=N1, N2=N2,
                deltaz=deltaz)
        c.x3 = c._x3_CST
        return c

    @classmethod
    def cylindrical(cls, D, chord=1, color='k', configuration=0):
        c = cls(D, color=color, chord=chord, n=1, configuration=configuration)
        c.x3 = c._x3_cylindrical
        return c

    def _x3_poly(self, x1, diff=None, D=None):
        """ z2 (checked)"""
        if D is None:
            D = self.D
        if diff is None:
            return(D[5]*x1**5 + D[4]*x1**4 + D[3]*x1**3 + D[2]*x1**2 + D[1]*x1 + D[0])
        elif diff == 'x1':
            return(5*D[5]*x1**4 + 4*D[4]*x1**3 + 3*D[3]*x1**2 + 2*D[2]*x1 + D[1])
        elif diff == 'x11':
            return(20*D[5]*x1**3 + 12*D[4]*x1**2 + 6*D[3]*x1 + 2*D[2])
        elif diff == 'x111':
            return(60*D[5]*x1**2 + 24*D[4]*x1 + 6*D[3])
        elif diff == 'theta1':
            return self.x3(x1, 'x1')*self.x1(x1, 'theta1')
        elif diff == 'theta11':
            # print('x3/x11: ', self.x3(x1, 'x11'))
            # print('x1/t1: ', self.x1(x1, 'theta1'))
            # print('x3/x1: ', self.x3(x1, 'x1'))
            # print('x1/t1: ', self.x1(x1, 'theta11'))
            return self.x3(x1, 'x11')*self.x1(x1, 'theta1')**2 + \
                self.x3(x1, 'x1')*self.x1(x1, 'theta11')
        elif diff == 'theta3':
            return(self.a[2, :, 2])
            # return(np.ones(len(x1)))
        elif diff == 'theta31' or diff == 'theta13':
            return(self.x1(x1, 'theta11'))
            # return(np.zeros(len(x1)))
        else:
            return(np.zeros(len(x1)))

    def _x3_CST(self, x1, diff=None):
        A = self.D[:-1]
        psi = x1 / self.chord
        if diff is None:
            return(CST(x1, self.chord, deltasz=self.deltaz, Au=A,
                       N1=self.N1, N2=self.N2))
        elif diff == 'x1':
            return(dxi_u(psi, A, self.deltaz/self.chord,
                         N1=self.N1, N2=self.N2))
        elif diff == 'x11':
            return(ddxi_u(psi, A, N1=self.N1, N2=self.N2))
        elif diff == 'theta1':
            return self.x3(x1, 'x1')*self.x1(x1, 'theta1')
        elif diff == 'theta11':
            return self.x3(x1, 'x11')*self.x1(x1, 'theta1')**2 + \
                self.x3(x1, 'x1')*self.x1(x1, 'theta11')
        elif diff == 'theta3':
            return(np.zeros(len(x1)))

    def _x3_cylindrical(self, x1, diff=None, R=None):
        """ z2 (checked)"""
        if R is None:
            R = self.D[0]
        # print(x1, R)
        if diff is None:
            return(R - np.sqrt(R**2 - x1**2))
        elif diff == 'x1':
            return(x1/np.sqrt(R**2 - x1**2))
        elif diff == 'x11':
            return(R**2/np.sqrt(R**2 - x1**2)**3)
        elif diff == 'theta3':
            return(np.zeros(len(x1)))

    def r(self, x1=None, diff=None):
        if x1 is None:
            x1 = self.x1_grid
        else:
            if type(x1) == float:
                x1 = np.array([x1])
        if diff == 'theta1':
            output = np.array([self.x1(x1, 'x1'),
                               self.x3(x1, 'x1')]).T
            output = np.einsum('ij,i->ij', output, self.x1(x1, 'theta1'))
        else:
            output = np.array([self.x1(x1, diff),
                               self.x3(x1, diff)]).T
            self.position = output
        return (output)

    def basis(self, x1=None, diff=None):
        if x1 is None:
            x1 = self.x1_grid

        if diff is None:
            self.a = np.zeros([3, len(x1), 3])
            self.a[0, :, :] = np.array([self.x1(x1, 'x1')*self.x1(x1, 'theta1'),
                                        [0]*len(x1),
                                        self.x3(x1, 'x1')*self.x1(x1, 'theta1')]).T
            # self.a[0,:,:] = np.array([[1]*len(x1),
            #                          [0]*len(x1),
            #                           [0]*len(x1)]).T
            self.a[1, :, :] = np.array([[0, 1, 0], ]*len(x1))
            self.a[2, :, :] = np.cross(self.a[0, :, :], self.a[1, :, :])

        elif diff == 'theta':
            # Most components are null
            self.da = np.zeros([3, 3, len(x1), 3])
            #  a1 diff theta1
            for i in range(3):
                for j in range(3):

                    # cross terms are null for 2D case
                    x11_1 = np.einsum('ij,i->ij', self.r(x1, 'x11'), self.x1(x1,
                                                                             'theta%d' % (i+1))*self.x1(x1, 'theta%d' % (j+1)))
                    x11_2 = np.einsum('ij,i->ij', self.r(x1, 'x1'),
                                      self.x1(x1, 'theta%d%d' % (i+1, j+1)))
                    x22_1 = np.einsum('ij,i->ij', self.r(x1, 'x22'), self.x2(x1,
                                                                             'theta%d' % (i+1))*self.x1(x1, 'theta%d' % (j+1)))
                    x22_2 = np.einsum('ij,i->ij', self.r(x1, 'x2'),
                                      self.x2(x1, 'theta%d%d' % (i+1, j+1)))
                    self.da[i, j, :, :] = x11_1+x11_2 + x22_1+x22_2

    def christoffel(self, i, j, k, order=1):
        if order == 1:
            gik_j = self.dA[i, k, j]
            gjk_i = self.dA[j, k, i]
            gij_k = self.dA[i, j, k]
            # print('dA', i,k,j, gik_j)
            # print('dA', j,k,i, gjk_i)
            # print('dA', i,j,k, gij_k)
            return .5*(gik_j + gjk_i - gij_k)
        elif order == 2:
            raise NotImplementedError

    def metric_tensor(self, diff=None):

        if diff is None:
            self.A = np.zeros([3, 3, len(self.x1_grid)])
            for i in range(3):
                for j in range(3):
                    self.A[i, j] = np.einsum('ij,ij->i', self.a[i, :], self.a[j, :])

        elif diff == 'theta':
            self.dA = np.zeros([3, 3, 3, len(self.x1_grid)])
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        self.dA[i, j, k] = np.einsum('ij,ij->i', self.da[i, k, :],
                                                     self.a[j, :]) + \
                            np.einsum('ij,ij->i', self.a[i, :],
                                      self.da[j, k, :])
                        # if self.dA[i,j,k,0] !=0:
                        #     print(i,j,k)
                        #     print(self.da[i,k,:])
                        #     print(self.a[j,:])
                        #     print(self.a[i,:])
                        #     print(self.da[j,k,:])

    def curvature_tensor(self):
        self.B = np.zeros([2, 2, len(self.x1_grid)])
        for alpha in range(2):
            for beta in range(2):
                # print(alpha, beta)
                # print(self.christoffel(alpha, beta, 2))
                self.B[alpha, beta] = self.christoffel(alpha, beta, 2)

    def arclength(self, chord=None):
        def integrand(x1):
            # dr = self.r(x1, 'x1')
            # if np.isnan(np.sqrt(np.inner(dr, dr))):
            #     return(100)
            # else:
            #     return np.sqrt(np.inner(dr, dr)[0, 0])
            dr = self.x3(np.array([x1]), 'x1')
            if np.isnan(dr):
                # print('NaN', x1)
                if x1 == 0:
                    dr = self.x3(np.array([1e-6]), 'x1')
                else:
                    dr = self.x3(np.array([x1-1e-6]), 'x1')
            return np.sqrt(1 + dr**2)
        if chord is None:
            chord = self.chord
        return integrate.quad(integrand, 0, chord, limit=500)

    def calculate_x1(self, length_target, bounds=None, output=False):
        def f(c_c):
            length_current, err = self.arclength(c_c[0])
            return abs(target - length_current)
        x0 = 0
        x1 = []

        for target in length_target:
            # print(target, x0)
            x1.append(optimize.fsolve(f, x0)[0])
            x0 = x1[-1]
        if output:
            return np.array(x1)
        else:
            self.x1_grid = np.array(x1)

    def plot(self, basis=False, r=None, label=None, linestyle='-', color=None, scatter=False, zorder=0, marker='.'):
        if r is None:
            r = self.r(self.x1_grid)

        if color is None:
            color = self.color

        if scatter:
            if label is None:
                plt.scatter(r[:, 0], r[:, 1], c=color, zorder=2, marker=marker)
            else:
                plt.scatter(r[:, 0], r[:, 1], c=color, label=label,
                            zorder=2, edgecolors='k', marker=marker)
        else:
            if label is None:
                plt.plot(r[:, 0], r[:, 1], color, linestyle=linestyle, lw=3,
                         zorder=zorder)
            else:
                plt.plot(r[:, 0], r[:, 1], color, linestyle=linestyle, lw=3,
                         label=label, zorder=zorder)
        if basis:
            plt.quiver(r[:, 0], r[:, 2],
                       self.a[0, :, 0], self.a[0, :, 2],
                       angles='xy', color=color, scale_units='xy')
            plt.quiver(r[:, 0], r[:, 2],
                       self.a[2, :, 0], self.a[2, :, 2],
                       angles='xy', color=color, scale_units='xy')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')

    def radius_curvature(self, x):
        self.rho = self.x3(x, diff='x11')/(1+(self.x3(x, diff='x1'))**2)**(3/2)

    def internal_variables(self, target_length):
        # At first we utilize the non-dimensional trailing edge thickness
        self.deltaz = self.D[-1]
        self.chord = 1
        nondimensional_length, err = self.arclength(chord=1.)
        self.chord = target_length/nondimensional_length
        self.deltaz = self.deltaz*self.chord
