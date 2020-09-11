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
        if hasattr(self, 'name') and self.name == 'pCST':
            self._D = values
            # Update internal variables if class already populated
            if hasattr(self, 'cst'):
                self._pCST_update()

        else:
            if self.n == 1 and (isinstance(values, float) or isinstance(values, np.float64)):
                self._D = values
            else:
                if len(values) != self.n:
                    self._D = np.zeros(self.n)
                    self._D[:len(values)] = values
                else:
                    self._D = values

    @property
    def x1_grid(self):
        return self._x1_grid

    @x1_grid.setter
    def x1_grid(self, values):
        if type(values) == list:
            values = np.array(values)

        if self.name == 'pCST':
            c_min = 0
            for i in range(self.p):
                c_max = c_min + self.cst[i].chord
                indexes = np.where((values >= (c_min-self.tol)) & (values <= (c_max+self.tol)))
                self.cst[i].x1_grid = values[indexes]
                c_min = c_max
        self._x1_grid = values

    def _x1(self, x1, diff=None):
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

    def _x1_pCST(self, x1, diff=None):
        output = np.zeros(len(x1))
        c_min = 0
        for i in range(self.p):
            c_max = c_min + self.cst[i].chord
            indexes = np.where((x1 >= (c_min-self.tol)) & (x1 <= (c_max+self.tol)))
            output[indexes] = self.cst[i].x1(x1[indexes], diff=diff)
            c_min = c_max
        return output

    def x2(self, x1, x2_value=0, diff=None):
        if diff is None:
            return np.zeros(len(x1))
        if diff is 'x2' or diff is 'theta2':
            return np.ones(len(x1))

    @classmethod
    def polynomial(cls, D, chord=1, color='b', n=6, tol=1e-6):
        c = cls(D, chord=chord, color=color, n=n, tol=tol, offset_x=0)
        c.x1 = c._x1
        c.x3 = c._x3_poly
        c.name = 'polynomial'
        return c

    @classmethod
    def CST(cls, D, chord, color, N1=.5, N2=1.0, deltaz=0, tol=1e-6, offset=0,
            deltazLE=0, offset_x=0):
        c = cls(D, chord=chord, color=color, n=len(D), N1=N1, N2=N2,
                deltaz=deltaz, tol=tol, offset=offset, deltazLE=deltazLE,
                offset_x=offset_x)
        c.x1 = c._x1
        c.x3 = c._x3_CST
        c.name = 'CST'
        c.zetaT = c.deltaz/c.chord
        c.zetaL = c.deltazLE/c.chord
        c.length = c.arclength()[0]
        return c

    @classmethod
    def pCST(cls, D, chord=np.array([.5, .5]), color='b', N1=[1, 1.], N2=[1.0, 1.0],
             tol=1e-6, offset=0, offset_x=0):
        c = cls(D, chord=chord, color=color, n=len(D), N1=N1, N2=N2,
                tol=tol, offset=offset, offset_x=offset_x)
        c.x1 = c._x1_pCST
        c.x3 = c._x3_pCST
        c.name = 'pCST'
        c.p = len(chord)
        n = int((len(D)-2)/c.p)

        c.cst = []
        c.zetaL = []
        c.zetaT = []
        c.A0 = []
        offset_s = 0
        for i in range(c.p):
            j = i - 1
            # From shape coefficients 1 to n
            Ai = D[1+i*n:1+(i+1)*n]
            Aj = D[1+j*n:1+(j+1)*n]
            if i == 0:
                offset_x = 0
                c.A0.append(D[0])
                c.zetaT.append(D[-1])
                c.zetaL.append(0)
            else:
                if N1[i] == 1. and N2[i] == 1.:
                    offset_x = chord[j]
                    ddj = n*Aj[-2] - (N1[j]+n)*Aj[-1]
                    c.A0.append((-chord[i]/chord[j]*ddj+Ai[0]*n)/(n+1))
                    c.zetaL.append(chord[j]/chord[i]*c.zetaT[j])
                    c.zetaT.append(-Aj[-1] + c.zetaT[j] - c.A0[-1] + c.zetaL[i] - c.zetaL[j])
                else:
                    raise(NotImplementedError)
            Di = [c.A0[i]] + list(Ai) + [c.zetaT[i]]
            c.cst.append(CoordinateSystem.CST(Di, chord[i],
                                              color[i], N1=N1[i], N2=N2[i],
                                              deltaz=c.zetaT[-1]*chord[i],
                                              deltazLE=c.zetaL[-1]*chord[i],
                                              offset_x=offset_x, offset=offset))
            c.cst[i].offset_s = offset_s
            offset_s += c.cst[i].length
        c.total_chord = sum([c.cst[i].chord for i in range(c.p)])
        c.total_length = sum([c.cst[i].length for i in range(c.p)])
        return c

    def _pCST_update(self):
        n = int((len(self.D)-2)/self.p)
        offset_s = 0
        for i in range(self.p):
            j = i - 1
            # From shape coefficients 1 to n
            Ai = self.D[1+i*n:1+(i+1)*n]
            Aj = self.D[1+j*n:1+(j+1)*n]
            error = 999
            while error > 1e-4:
                chord0 = self.cst[i].chord
                if i == 0:
                    offset_x = 0
                    self.A0[i] = self.D[0]
                    self.zetaT[i] = self.D[-1]
                    self.zetaL[i] = 0
                else:
                    if self.N1[i] == 1. and self.N2[i] == 1.:
                        offset_x = self.cst[j].chord
                        ddj = n*Aj[-2] - (self.N1[j]+n)*Aj[-1]
                        self.A0[i] = (-self.cst[i].chord/self.cst[j].chord*ddj+Ai[0]*n)/(n+1)
                        self.zetaL[i] = self.cst[j].chord/self.cst[i].chord*self.zetaT[j]
                        self.zetaT[i] = -Aj[-1] + self.zetaT[j] - self.A0[i] + \
                            self.zetaL[i] - self.zetaL[j]
                    else:
                        raise(NotImplementedError)
                Di = [self.A0[i]] + list(Ai) + [self.zetaT[i]]
                self.cst[i].D = Di
                self.cst[i].zetaT = self.zetaT[i]
                self.cst[i].zetaL = self.zetaL[i]
                self.cst[i].offset_x = offset_x
                self.cst[i].internal_variables(self.cst[i].length)
                if i == 0:
                    error = 0
                else:
                    error = abs(chord0 - self.cst[i].chord)
            self.cst[i].offset_s = offset_s
            offset_s += self.cst[i].length
        self.total_chord = sum([self.cst[i].chord for i in range(self.p)])

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

    def _x3_CST(self, x1, diff=None, offset=True):
        A = self.D[:-1]
        if offset:
            x1 = x1 - self.offset_x
        psi = x1 / self.chord
        if diff is None:
            return(self.offset + CST(x1, self.chord, deltasz=self.zetaT*self.chord, Au=A,
                                     N1=self.N1, N2=self.N2, deltasLE=self.zetaL*self.chord))
        elif diff == 'x1':
            d = dxi_u(psi, A, self.zetaT, N1=self.N1, N2=self.N2, zetaL=self.zetaL)
            if x1[0] < 1e-6 and self.N1 == 1:
                d[0] = +A[0] + self.zetaT - self.zetaL
            if abs(x1[-1] - self.chord) < 1e-6 and self.N2 == 1:
                d[-1] = -A[-1] + self.zetaT - self.zetaL
            return d
        elif diff == 'x11':
            dd = (1/self.chord)*ddxi_u(psi, A, N1=self.N1, N2=self.N2)
            if abs(x1[0]) < 1e-6 and self.N1 == 1 and self.N2 == 1:
                n = len(A) - 1
                dd[0] = (1/self.chord)*(-2*(n+1)*A[0] + 2*A[1]*n)
            if abs(x1[-1] - self.chord) < 1e-6 and self.N2 == 1:
                n = len(A) - 1
                dd[-1] = (1/self.chord)*(2*n*A[-2] - 2*(self.N1+n)*A[-1])
            return dd
        elif diff == 'theta1':
            return self.x3(x1, 'x1')*self.x1(x1, 'theta1')
        elif diff == 'theta11':
            return self.x3(x1, 'x11')*self.x1(x1, 'theta1')**2 + \
                self.x3(x1, 'x1')*self.x1(x1, 'theta11')
        elif diff == 'theta3':
            return(np.zeros(len(x1)))

    def _x3_pCST(self, x1, diff=None):
        output = np.zeros(len(x1))
        c_min = 0
        for i in range(self.p):
            c_max = c_min + self.cst[i].chord
            indexes = np.where((x1 >= (c_min-self.tol)) & (x1 <= (c_max+self.tol)))
            output[indexes] = self.cst[i].x3(x1[indexes], diff=diff)
            c_min = c_max
        return output

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

    def arclength(self, chord=None, origin=None):
        def integrand(x1):
            dr = self.x3(np.array([x1]), 'x1', offset=False)
            return np.sqrt(1 + dr**2)

        if chord is None:
            chord = self.chord
        if origin is None:
            origin = self.origin
        if chord == 0:
            return [0, 0]
        return integrate.quad(integrand, origin, chord, limit=500)[0]

    def arclength_index(self, index):
        x1 = self.x1_grid[index]
        dr = self.x3(np.array([x1]), 'x1')
        if np.isnan(dr):
            if x1 == 0:
                dr = self.x3(np.array([self.tol]), 'x1')
            else:
                dr = self.x3(np.array([x1-self.tol]), 'x1')
        self.darc[index] = np.sqrt(1 + dr[0]**2)
        return integrate.trapz(self.darc[:index+1], self.x1_grid[:index+1]-self.offset_x)

    def improper_arclength_index(self, index):
        x1 = self.x1_grid[index]
        self.darc[index] = self.bounded_dr(x1)
        bounded_integral = integrate.trapz(self.darc[:index+1], self.x1_grid[:index+1])
        return bounded_integral + self.unbounded_integral(x1, 0)

    def arclength_chord(self):
        self.darc = np.ones(len(self.x1_grid))
        dr = self.x3(np.array([1e-7]), 'x1', offset=False)
        self.darc[0] = np.sqrt(1 + dr[0]**2)
        for index in range(1, len(self.x1_grid)-1):
            x1 = self.x1_grid[index]
            dr = self.x3(np.array([x1]), 'x1', offset=False)
            self.darc[index] = np.sqrt(1 + dr[0]**2)
        self.darc[-1] = np.sqrt(1 + (-self.D[-2]+self.D[-1]-self.zetaL)**2)
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
            dr = self.x3(np.array([x]), 'x1', offset=False)
            return np.sqrt(1 + dr[0]**2) - np.sqrt(1 + A/x)

    def unbounded_integral(self, end, start=0):
        def indefinite_integral(x):
            return np.sqrt(x*(A + x)) + A*np.log(np.sqrt(A+x) + np.sqrt(x))
        A = self.N1**2*self.D[0]**2
        return indefinite_integral(end) - indefinite_integral(start)

    def calculate_x1(self, length_target, bounds=None, output=False, origin=0, length_rigid=0):
        def f(c_c):
            length_current, err = self.arclength(c_c[0])
            # print(c_c, length_current)
            return abs(target - length_current)

        def f_index(x):
            # Penalize in case x goes negative
            if x < 0:
                return 100
            else:
                self.x1_grid[index] = x + self.offset_x
                # if self.name == 'CST':
                #     length_current = length_rigid + self.improper_arclength_index(index)
                # else:
                length_current = length_rigid + self.arclength_index(index)
                return abs(target - length_current)

        if self.name == 'pCST':
            self.s = length_target
            s_min = 0
            for i in range(self.p):
                s_max = s_min + self.cst[i].length
                indexes = np.where((length_target >= (s_min-self.tol)) &
                                   (length_target <= (s_max+self.tol)))
                self.cst[i].s = length_target[indexes]
                self.cst[i].calculate_x1(length_target[indexes] - self.cst[i].offset_s)
                s_min = s_max
            x1_grid = []
            for i in range(self.p):
                x1_grid += list(self.cst[i].x1_grid)
            self.x1_grid = x1_grid  # [self.cst[i].x1_grid for i in range(self.p)]
        else:
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
                        dr = self.x3(np.array([origin + self.offset_x]), 'x1')
                        self.darc[index] = np.sqrt(1 + dr[0]**2)
                        self.x1_grid[index] = origin + self.offset_x
                    else:
                        target = length_target[index]
                        self.x1_grid[index] = optimize.fsolve(f_index, target)[0] + self.offset_x
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
            s_epsilon = self.arclength(np.array([origin]))
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
        if self.name == 'pCST':
            for i in range(self.p):
                self.cst[i].plot(basis, label=label[i], scatter=scatter)
        else:
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
            if self.name == 'CST':
                if x[0] == 0:
                    if self.D[0] == 0:
                        rho[0] = 0
                    else:
                        rho[0] = -2/(self.D[0]**2)/self.chord
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
        # print('target', target_length, origin, self.D, self.N1,
        #       self.N2, target_length/nondimensional_length)
        self.chord = target_length/nondimensional_length
        self.deltaz = self.zetaT*self.chord
        self.deltazLE = self.zetaL*self.chord

    def calculate_angles(self):
        self.cos = self.x1(self.x1_grid, 'theta1')
        self.sin = self.x3(self.x1_grid, 'theta1')

        if self.x1_grid[0] == 0 and self.N1 == 1:
            self.cos[0] = 1
            self.sin[0] = 0
