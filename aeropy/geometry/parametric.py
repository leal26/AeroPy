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
    def polynomial(cls, D, chord=1, color='b', n=6, tol=1e-6):
        c = cls(D, chord=chord, color=color, n=n, tol=tol)
        c.x3 = c._x3_poly
        c.name = 'polynomial'
        return c

    @classmethod
    def CST(cls, D, chord, color, N1=.5, N2=1.0, deltaz=0, tol=1e-6):
        c = cls(D, chord=chord, color=color, n=len(D), N1=N1, N2=N2,
                deltaz=deltaz, tol=tol)
        c.x3 = c._x3_CST
        c.name = 'CST'
        c.zetaT = c.deltaz/c.chord
        return c

    @classmethod
    def cylindrical(cls, D, chord=1, color='k', configuration=0):
        c = cls(D, color=color, chord=chord, n=1, configuration=configuration)
        c.x3 = c._x3_cylindrical
        c.name = 'cylindrical'
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
            d = dxi_u(psi, A, self.zetaT, N1=self.N1, N2=self.N2)
            if abs(x1[-1] - self.chord) < 1e-5:
                d[-1] = -A[-1] + self.zetaT
            return d
        elif diff == 'x11':
            return((1/self.chord)*ddxi_u(psi, A, N1=self.N1, N2=self.N2))
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

    def curvature_tensor(self):
        self.B = np.zeros([2, 2, len(self.x1_grid)])
        for alpha in range(2):
            for beta in range(2):
                self.B[alpha, beta] = self.christoffel(alpha, beta, 2)

    def arclength(self, chord=None, origin=0):
        def integrand(x1):
            dr = self.x3(np.array([x1]), 'x1')
            return np.sqrt(1 + dr**2)
        if chord is None:
            chord = self.chord
        return integrate.quad(integrand, origin, chord, limit=500)

    def arclength_index(self, index):
        x1 = self.x1_grid[index]
        dr = self.x3(np.array([x1]), 'x1')
        if np.isnan(dr):
            if x1 == 0:
                dr = self.x3(np.array([self.tol]), 'x1')
            else:
                dr = self.x3(np.array([x1-self.tol]), 'x1')
        self.darc[index] = np.sqrt(1 + dr[0]**2)
        return integrate.trapz(self.darc[:index+1], self.x1_grid[:index+1])

    def improper_arclength_index(self, index):
        x1 = self.x1_grid[index]
        self.darc[index] = self.bounded_dr(x1)
        bounded_integral = integrate.trapz(self.darc[:index+1], self.x1_grid[:index+1])
        return bounded_integral + self.unbounded_integral(x1, 0)

    def arclength_chord(self):
        self.darc = np.ones(len(self.x1_grid))
        dr = self.x3(np.array([1e-7]), 'x1')
        self.darc[0] = np.sqrt(1 + dr[0]**2)
        for index in range(1, len(self.x1_grid)-1):
            x1 = self.x1_grid[index]
            dr = self.x3(np.array([x1]), 'x1')
            self.darc[index] = np.sqrt(1 + dr[0]**2)
        self.darc[-1] = np.sqrt(1 + (-self.D[-2]+self.D[-1])**2)
        return integrate.trapz(self.darc, self.x1_grid)

    def improper_arclength_chord(self):
        self.darc = np.zeros(len(self.x1_grid))
        for index in range(1, len(self.x1_grid)):
            x1 = self.x1_grid[index]
            self.darc[index] = self.bounded_dr(x1)
        bounded_integral = integrate.trapz(self.darc, self.x1_grid)
        return bounded_integral + self.unbounded_integral(self.chord)

    def bounded_dr(self, x):
        A = self.N1**2*self.D[0]**2
        n = self.n - 2
        if x == 0:
            return self.D[0]**2*(-n-1) + 1.5*n*self.D[0]*self.D[1]
        elif x == 1:
            # - np.sqrt(1 + B/np.sqrt(x))
            return np.sqrt(1 + (-self.D[-2]+self.D[-1])**2) - np.sqrt(1 + A/x)
        else:
            dr = self.x3(np.array([x]), 'x1')
            return np.sqrt(1 + dr[0]**2) - np.sqrt(1 + A/x)

    def unbounded_integral(self, end, start=0):
        def indefinite_integral(x):
            return np.sqrt(x*(A + x)) + A*np.log(np.sqrt(A+x) + np.sqrt(x))
        A = self.N1**2*self.D[0]**2
        return indefinite_integral(end) - indefinite_integral(start)

    def calculate_x1(self, length_target, bounds=None, output=False, origin=0,
                     length_rigid=0):
        def f(c_c):
            length_current, err = self.arclength(c_c[0])
            return abs(target - length_current)

        def f_index(dx):
            if dx[0] < 0:
                return 100
            else:
                self.x1_grid[index] = self.x1_grid[index-1] + dx[0]

                if self.name == 'CST':
                    length_current = length_rigid + self.improper_arclength_index(index)
                else:
                    length_current = length_rigid + self.arclength_index(index)
                return target - length_current

        def fprime_index(x):
            dr = self.x3(x, 'x1')
            return np.array([np.sqrt(1 + dr[0]**2)])

        if len(length_target) == 1:
            target = length_target[0]
            x1 = [optimize.fsolve(f, 0, fprime=fprime_index)[0]]
        else:
            x1 = [origin]
            if hasattr(self, 'x1_grid'):
                if np.isnan(self.x1_grid).any():
                    prev_values = False
                else:
                    prev_values = True
            else:
                prev_values = False
            if not prev_values:
                self.x1_grid = np.zeros(len(length_target))
                self.darc = np.ones(len(length_target))
                self.x1_grid[0] = origin
                if not self.name == 'CST' and origin != 0:
                    dr = self.x3(np.array([origin]), 'x1')
                    self.darc[0] = np.sqrt(1 + dr**2)
                elif self.name == 'CST' and origin != 0:
                    raise(NotImplementedError)
            for index in range(1, len(length_target)):
                target = length_target[index]
                if prev_values:
                    x0 = self.x1_grid[index] - self.x1_grid[index-1]
                else:
                    x0 = length_target[index] - length_target[index-1]
                dx = optimize.fsolve(f_index, x0, fprime=fprime_index)[0]
                x1.append(self.x1_grid[index-1] + dx)
        if output:
            return np.array(x1)
        else:
            self.x1_grid = np.array(x1)

    def calculate_x1(self, length_target, bounds=None, output=False, origin=0, length_rigid=0):
        def f(c_c):
            length_current, err = self.arclength(c_c[0])
            return abs(target - length_current)

        def f_index(x):
            # Penalize in case x goes negative
            if x < 0:
                return 100
            else:
                self.x1_grid[index] = x
                if self.name == 'CST':
                    length_current = length_rigid + self.improper_arclength_index(index)
                else:
                    length_current = length_rigid + self.arclength_index(index)
                return abs(target - length_current)

    #     def fprime_index(x):
    #         dr = self.x3(x, 'x1')
    #         return np.array([np.sqrt(1 + dr[0]**2)])
        x1 = []

        if len(length_target) == 1:
            target = length_target[0]
            x1.append(optimize.fsolve(f, origin)[0])
        else:
            self.x1_grid = np.zeros(len(length_target))
            self.darc = np.zeros(len(length_target))
            for index in range(len(length_target)):
                if index == 0 and origin != 0:
                    dr = self.x3(np.array([origin]), 'x1')
                    self.darc[index] = np.sqrt(1 + dr[0]**2)
                    self.x1_grid[index] = origin
                else:
                    target = length_target[index]
                    self.x1_grid[index] = optimize.fsolve(f_index, target)[0]
        if output:
            return np.array(x1)

    def calculate_s(self, N, target_length=None, density='gradient', origin=0):
        def integrand(s):
            rho = self.radius_curvature(np.array([s]), output_only=True)[0]
            return abs(rho)

        def f(dx):
            if dx[0] >= 0:
                self.x1_grid[i] = self.x1_grid[i-1] + dx[0]
                self.rho[i] = integrand(self.x1_grid[i])
                partial = integrate.quad(integrand,
                                         self.x1_grid[i-1], self.x1_grid[i-1]+dx[0], limit=500)[0]
                self.partial[i] = partial
                return abs(partial - total/(N-1))
            else:
                return 100

        def fprime(dx):
            x = self.x1_grid[i-1]+dx[0]
            output = integrand(x)
            return [output]

        if density == 'gradient':
            s_epsilon = self.arclength(np.array([origin]))[0]
            return np.linspace(s_epsilon, target_length, N)
        elif density == 'curvature':
            total = integrate.quad(integrand, origin, self.chord, limit=500)[0]

            self.x1_grid = np.zeros(N)
            self.darc = np.zeros(N)
            self.partial = np.zeros(N)
            self.rho = np.zeros(N)
            self.rho[0] = abs(self.radius_curvature(np.array([0]), output_only=True)[0])
            s_list = [0]

            for i in range(1, N):
                dx = optimize.fsolve(f, 0, fprime=fprime)[0]
                self.x1_grid[i] = self.x1_grid[i-1] + dx
                s_list.append(self.improper_arclength_index(i))
                # BREAK
            return np.array(s_list)

    def plot(self, basis=False, r=None, label=None, linestyle='-', color=None, scatter=False, zorder=0, marker='.'):
        if r is None:
            r = self.r(self.x1_grid)
        print('x_p', r[:, 0])
        print('y_p', r[:, 1])
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

    def radius_curvature(self, x, output_only=False, parametric=False):
        if parametric:
            if self.name == 'CST' and x[0] == 0:
                x[0] = 1e-7
            dx = self.x1(x, diff='theta1')
            dy = self.x3(x, diff='theta1')
            ddx = self.x1(x, diff='theta11')
            ddy = self.x3(x, diff='theta11')
            rho = (dx*ddy-dy*ddx)/(dx**2 + dy**2)**(1.5)
        else:
            rho = self.x3(x, diff='x11')/(1+(self.x3(x, diff='x1'))**2)**(3/2)
            if self.name == 'CST' and x[0] == 0:
                if self.D[0] == 0:
                    rho[0] = 0
                else:
                    rho[0] = -2/(self.D[0]**2)/self.chord
                # x_i = np.array([1e-7])
                # r = self.x3(x_i, diff='x11')/(1+(self.x3(x_i, diff='x1'))**2)**(3/2)
                # rho[0] = r
        if output_only:
            return rho
        else:
            self.rho = rho

    def internal_variables(self, target_length, origin=0):
        # At first we utilize the non-dimensional trailing edge thickness
        self.zetaT = self.D[-1]
        origin = origin/self.chord
        self.chord = 1

        nondimensional_length, err = self.arclength(chord=1., origin=origin/self.chord)

        self.chord = target_length/nondimensional_length
        self.deltaz = self.zetaT*self.chord

    def calculate_angles(self):
        self.cos = self.x1(self.x1_grid, 'theta1')
        self.sin = self.x3(self.x1_grid, 'theta1')

        if self.x1_grid[0] == 0 and self.N1 == 1:
            self.cos[0] = 1
            self.sin[0] = 0
