import aeropy
import math
import copy

from findiff import FinDiff
import numpy as np
from scipy.integrate import quad, trapz
from scipy.optimize import minimize, curve_fit
import matplotlib.pyplot as plt


class euler_bernoulle_curvilinear():
    def __init__(self, geometry_parent, geometry_child, properties,
                 load, load_type, theta1):
        self.properties = properties
        self.load = load
        self.load_type = load_type

        # Defining geometries
        self.g_p = geometry_parent
        self.g_c = geometry_child

        # updating s CoordinateSystem
        self.theta1 = theta1
        self.g_p.calculate_x1(self.theta1)
        self.g_c.calculate_x1(self.theta1)
        self.arc_length = self.g_p.arclength()[0]

    def update_chord(self, length_target=None, bounds=None):
        def f(c_c):
            length_current, err = self.g_c.arclength(c_c)
            return abs(length_target - length_current)
        if length_target is None:
            length_target = self.arc_length
        if bounds is None:
            self.g_c.chord = minimize(f, self.g_c.chord).x[0]
        else:
            self.g_c.chord = minimize(f, self.g_c.chord,
                                      method='L-BFGS-B',
                                      bounds=bounds).x[0]
        # In case the calculated chord is really close to the original
        if abs(self.g_p.chord - self.g_c.chord) < 1e-7:
            self.g_c.chord = self.g_p.chord

    def analytical_solutions(self):
        if self.load_type == 'concentrated':
            bp = self.properties
            self.g_c.D = self.load/(6*bp.young*bp.inertia) * \
                np.array([0, 0, 3, -1, 0])
        elif self.load_type == 'distributed':
            bp = self.properties
            c2 = 6*(bp.length**2)*self.load/(24*bp.young*bp.inertia)
            c3 = -4*(bp.length)*self.load/(24*bp.young*bp.inertia)
            c4 = (1)*self.load/(24*bp.young*bp.inertia)
            self.g_c.D = np.array([0, 0, c2, c3, c4])

    def bending_strain(self):
        self.B = self.g_c.x3(self.g_c.x1_grid, 'x11') - \
            self.g_p.x3(self.g_p.x1_grid, 'x11')

    def free_energy(self):
        bp = self.properties
        self.phi = bp.young*bp.inertia/2*self.B**2

    def strain_energy(self):
        self.U = np.trapz(self.phi, self.theta1)

    def work(self):
        u = self.g_c.x3(self.g_c.x1_grid) - self.g_p.x3(self.g_p.x1_grid)
        if self.load_type == 'concentrated':
            self.W = self.load*u[-1]
        elif self.load_type == 'distributed':
            self.W = np.trapz(np.multiply(self.load, u), self.g_c.x1_grid)

    def residual(self):
        self.R = self.U - self.W

    def minimum_potential(self, x0=[0, 0]):
        def to_optimize(x):
            self.g_c.D = [0, 0] + list(x)
            self.update_chord()
            self.g_c.calculate_x1(self.theta1)
            self.work()
            self.bending_strain()
            self.free_energy()
            self.strain_energy()
            self.residual()
            print(x, self.R)
            return self.R

        # With bounds
        bounds = np.array(((-0.1, 0.1),)*len(x0))

        res = minimize(to_optimize, x0)
        self.g_c.D = [0, 0] + list(res.x)
        self.work()
        self.free_energy()
        self.strain_energy()
        self.residual()
        return(res.x, res.fun)


class beam_chen():
    def __init__(self, geometry, properties, load, s, ignore_ends=False,
                 rotated=False, origin=0, g2=None, length_preserving=True,
                 calculate_resultant=False):
        self.g = copy.deepcopy(geometry)
        self.p = properties
        self.l = load
        if s is None:
            self.s = self.g.s
        else:
            self.s = s
        self.origin = origin
        self.ignore_ends = ignore_ends
        self.rotated = rotated
        self.length_preserving = length_preserving
        if self.ignore_ends:
            self.integral_ends()
        self.g_p = copy.deepcopy(geometry)
        # if length_preserving:
        #     self.g_p.internal_variables(self.g.length, origin=self.origin)
        self.g_p.calculate_x1(self.s, origin=origin, length_rigid=self.s[0])
        self.g_p.radius_curvature(self.g_p.x1_grid)

        if self.rotated:
            self.g_p.calculate_angles()
        self.calculate_resultant = calculate_resultant

    def calculate_M(self):
        self.M = np.zeros(len(self.s))
        self.repeat = False
        if self.g.name == 'pCST' and self.g.rigid_LE:
            for i in range(len(self.M)):
                j = i + len(self.g.cst[0].x1_grid)
                self.M[i] = self._M(self.x[j], self.s[i], self.y[j])
        else:
            for i in range(len(self.M)):
                self.M[i] = self._M(self.x[i], self.s[i], self.y[i])

    def _M(self, x, s, y=None):
        M_i = 0

        if self.l.concentrated_load is not None:
            if self.l.follower and self.rotated:
                for i in range(len(self.l.concentrated_s)):
                    index = np.where(self.s == self.l.concentrated_s[i])[0][0]
                    c = self.g_p.cos[index]*self.g.cos[index] + \
                        self.g_p.sin[index]*self.g.sin[index]
                    s = self.g_p.sin[index]*self.g.cos[index] - \
                        self.g_p.cos[index]*self.g.sin[index]
                    c = self.l.concentrated_direction[i][0]*self.g_p.cos[index] + \
                        self.l.concentrated_direction[i][1]*self.g_p.sin[index]
                    s = np.sqrt(1-c**2)

                    f2 = c*self.g.sin[index] - s*self.g.cos[index]
                    f1 = (c-self.g.sin[index]*f2)/self.g.cos[index]

                    M_i += self.l.concentrated_magnitude[i]*f1*(self.y[index]-y-y0)
                    M_i += self.l.concentrated_magnitude[i]*f2*(self.x[index]-x-x0)
            else:
                for i in range(len(self.l.concentrated_s)):
                    decimals = abs(int(math.floor(math.log10(abs(self.g.tol)))))
                    index = np.where(np.around(self.s, decimals=decimals) ==
                                     np.around(self.l.concentrated_s[i], decimals=decimals))[0][0]
                    if self.g.rigid_LE:
                        index += len(self.g.cst[0].x1_grid)
                    if s <= self.l.concentrated_s[i] + self.g.tol:
                        M_i -= self.l.concentrated_load[i][0]*(self.y[index]-y)
                        M_i += self.l.concentrated_load[i][1]*(self.x[index]-x)

        if self.l.distributed_load is not None:
            index = np.where(np.around(self.s, decimals=7) == np.around(s, decimals=7))[0][0]

            if self.ignore_ends:
                i_end = -1
            else:
                i_end = len(self.s)
            if not self.l.follower:
                M_i -= trapz(self.l.distributed_load(self.s[index:i_end])
                             * (self.x[index:i_end]-x-x0), self.s[index:i_end])
            else:
                w = self.l.distributed_load(self.s[index:i_end])
                M_x = w*self.g.cos[index:i_end]*(self.x[index:i_end]-x-x0)
                M_y = w*self.g.sin[index:i_end]*(self.y[index:i_end]-y-y0)
                M_i -= trapz(M_x + M_y, self.s[index:i_end])
        # Point torques
        for j in range(len(self.l.torque)):
            if np.around(s, decimals=7) <= self.l.torque_s[j]:
                if np.around(s, decimals=7) == self.l.torque_s[j]:
                    if self.repeat:
                        self.repeat = False
                    else:
                        self.repeat = True
                        M_i += self.l.torque[j]
                else:
                    M_i += self.l.torque[j]
        return M_i

    def calculate_G(self):
        self.G = np.zeros(len(self.x))
        for i in range(len(self.G)):
            self.G[i] = self._G(self.x[i], self.s[i])

    def _G(self, x, s):
        # deformed curvature radius
        c = 1/self.p.young/self.p.inertia
        index = np.where(self.x == x)[0][0]
        # return c*trapz(self.M[:index+1]*, self.s[:index+1])
        return c*trapz(self.M[:index+1], self.x[:index+1])

    def calculate_x(self):
        if not hasattr(self, 'x'):
            self.x = np.zeros(len(self.s))
        for i in range(len(self.s)):
            self.x[i] = self._x(self.s[i])

    def _x(self, s):
        def _to_minimize(l):
            index = np.where(self.s == s)[0][0]
            x = self.x[:index+1]
            x[-1] = l
            current_L = trapz(np.sqrt(1-self.G[:index+1]**2), x)
            return abs(s-current_L)
        return minimize(_to_minimize, s, method='Nelder-Mead',).x[0]

    def calculate_deflection(self):
        self.y = np.zeros(len(self.x))
        for i in range(len(self.x)):
            dydx = self.G[:i+1]/(1-self.G[:i+1]**2)
            y_i = trapz(dydx, self.x[:i+1])
            self.y[i] = y_i

    def calculate_residual(self):
        self.r = np.zeros(len(self.s))
        for i in range(len(self.M)):
            lhs = self.M[i]/self.p.young/self.p.inertia
            if self.g.rigid_LE:
                j = i + len(self.g.cst[0].x1_grid)
                rhs = self.g.rho[j] - self.g_p.rho[j]
            else:
                rhs = self.g.rho[i] - self.g_p.rho[i]
            self.r[i] = abs(lhs - rhs)
        if self.ignore_ends:
            self.R = abs(trapz(self.r[1:-1], self.s[1:-1]))
        else:
            self.R = abs(trapz(self.r[:], self.s[:]))
        if np.isnan(self.R):
            print('ignore_ends', self.ignore_ends)
            print('r', self.r)
            print('x', self.g.x1_grid)
            print('s', self.s)
            print('rho', self.g.rho)
            self.R = 100
            # BREAK

        # print('R: ', self.R)

    def iterative_solver(self):
        self.x = np.copy(self.s)
        x_before = np.copy(self.x)
        y_before = self.g.x3(self.x)

        error = 1000
        while error > 1e-8:
            self.calculate_M()
            self.calculate_G()
            self.calculate_x()
            self.calculate_M()
            self.calculate_deflection()
            error = np.linalg.norm((self.x-x_before)**2 + (self.y-y_before)**2)
            x_before = np.copy(self.x)
            y_before = np.copy(self.y)
            print(error)

    def parameterized_solver(self, format_input=None, x0=None, constraints=(),):
        def formatted_residual(A):
            A = format_input(A, self.g, self.g_p)
            return self._residual(A)

        bounds = []
        margin = 0.8
        for xi in x0:
            bounds.append([xi-margin, xi+margin])

        sol = minimize(formatted_residual, x0, method='SLSQP', bounds=bounds,
                       constraints=constraints)
        self.g.D = format_input(sol.x, self.g, self.g_p)
        if self.length_preserving and self.g.name != 'pCST':
            self.g.internal_variables(self.g.length, origin=self.origin)
        # self.g.calculate_x1(self.s)
        self.g.calculate_x1(self.s)
        self.x = self.g.x1_grid
        self.y = self.g.x3(self.x)
        print('sol', self.g.D, sol.fun)

    def _residual(self, A):
        self.g.D = A
        if self.length_preserving and self.g.name != 'pCST':
            self.g.internal_variables(self.g.length, origin=self.origin)

        self.g.calculate_x1(self.s)
        self.x = self.g.x1_grid

        self.y = self.g.x3(self.x)

        if self.l.follower:
            self.g.calculate_angles()
        self.g.radius_curvature(self.g.x1_grid)
        if self.calculate_resultant:
            self.calculate_resultants()
        self.calculate_M()
        self.calculate_residual()
        print('R', self.R)
        return self.R

    def integral_ends(self):
        # Correct point
        origin = np.array([self.origin])
        if self.g.name == 'pCST':
            tip = np.array([self.g.total_chord])
        else:
            tip = np.array([self.g.chord])
        if np.isnan(self.g.x3(origin, diff='x1')[0]) or \
                np.isnan(self.g.x3(origin, diff='x11')[0]):
            self.s = np.insert(self.s, 1, origin + self.g.tol)

        if np.isnan(self.g.x3(tip, diff='x1')[0]) or \
                np.isnan(self.g.x3(tip, diff='x11')[0]):
            self.s = np.insert(self.s, -1, self.s[-1] - self.g.tol)
        self.g.calculate_x1(self.s)

    def calculate_resultants(self):
        def M_f(s, *Rys):
            self.l.concentrated_load = self.l.external_load.copy()
            for Ry in Rys:
                self.l.concentrated_load += [[0, Ry], ]
            self.calculate_M()
            if self.ignore_ends:
                return self.M[1:-1]
            else:
                return self.M
        self.g.calculate_angles()

        M = (self.p.young*self.p.inertia)*(self.g.rho - self.g_p.rho)
        s = self.s[1:-1]
        self.l.concentrated_s = self.l.external_s.copy() + self.spars_s
        if self.ignore_ends:
            M = M[1:-1]
            s = s[1:-1]
        popt, pcov = curve_fit(M_f, s, M, p0=[0]*len(self.l.concentrated_s))

        self.l.concentrated_load = self.l.external_load.copy()
        for Ry in popt:
            self.Rx = 0  # -2.614  # sbeta*R
            self.Ry = Ry  # -3.634  # cbeta*R
            self.l.concentrated_load += [[self.Rx, self.Ry], ]


class airfoil():
    def __init__(self, g_upper, g_lower, properties_upper, properties_lower,
                 load_upper, load_lower, s_upper, s_lower, ignore_ends=False,
                 rotated=False, origin=0, chord=1, zetaT=0):
        self.bu = beam_chen(g_upper, properties_upper, load_upper, s_upper,
                            origin=origin, ignore_ends=ignore_ends,
                            rotated=rotated)
        self.bl = beam_chen(g_lower, properties_lower, load_lower, s_lower,
                            origin=origin, ignore_ends=ignore_ends,
                            rotated=rotated)
        self.chord = chord
        self.zetaT = zetaT

    def calculate_x(self, s_upper=None, s_lower=None):
        if s_upper is None and s_lower is None:
            self.bu.g.calculate_x1(self.bu.s, origin=self.bu.origin,
                                   length_rigid=self.bu.s[0])
            self.bl.g.calculate_x1(self.bl.s, origin=self.bl.origin,
                                   length_rigid=self.bl.s[0])
            self.bu.x = self.bu.g.x1_grid
            self.bl.x = self.bl.g.x1_grid
        else:
            raise(NotImplementedError)

    def parameterized_solver(self, format_input=None, x0=None):
        def formatted_residual(A):
            [Au, Al] = format_input(A, self.bu.g, self.bu.g_p, self.bl.g, self.bl.g_p)
            self.calculate_x()
            self.bl.g.radius_curvature(self.bl.g.x1_grid)
            self.bl.calculate_resultants()
            # self.bu.l.concentrated_load[0][0] = -self.bl.Rx
            # self.bu.l.concentrated_load[0][1] = -1 - self.bl.Ry
            R = self.bu._residual(Au) + self.bl._residual(Al)
            print('R', self.bu.R, self.bl.R, Au)
            return R

        sol = minimize(formatted_residual, x0, method='SLSQP')
        self.bu.g.D, self.bl.g.D = format_input(
            sol.x, self.bu.g, self.bu.g_p, self.bl.g, self.bl.g_p)
        # self.bu.g.internal_variables(self.bu.length, origin=self.bu.origin)
        # self.bl.g.internal_variables(self.bl.length, origin=self.bl.origin)
        self.calculate_x()

        self.bu.y = self.bu.g.x3(self.bu.x)
        self.bl.y = self.bl.g.x3(self.bl.x)
        print('sol', self.bu.g.D, self.bl.g.D)


class coupled_beams():
    def __init__(self, g_upper, g_lower, properties_upper, properties_lower,
                 load_upper, load_lower, s_upper, s_lower, ignore_ends=False,
                 rotated=False, origin=0, chord=1, zetaT=0, spars_s=None,
                 spars_k=None, coupled_forces=True):
        self.bu = beam_chen(g_upper, properties_upper, load_upper, s_upper,
                            origin=origin, ignore_ends=ignore_ends,
                            rotated=rotated)
        self.bl = beam_chen(g_lower, properties_lower, load_lower, s_lower,
                            origin=origin, ignore_ends=ignore_ends,
                            rotated=rotated)
        self.spars_s = spars_s
        self.spars_k = spars_k
        self.coupled_forces = coupled_forces

    def calculate_x(self, s_upper=None, s_lower=None):
        if s_upper is None and s_lower is None:
            self.bu.g.calculate_x1(self.bu.s, origin=self.bu.origin,
                                   length_rigid=self.bu.s[0])
            self.bl.g.calculate_x1(self.bl.s, origin=self.bl.origin,
                                   length_rigid=self.bl.s[0])
            self.bu.x = self.bu.g.x1_grid
            self.bl.x = self.bl.g.x1_grid
        else:
            raise(NotImplementedError)

    def parameterized_solver(self, format_input=None, x0=None, constraints=()):
        sol = minimize(self.formatted_residual, x0,  method='SLSQP', constraints=constraints, args=(format_input),
                       options={'iprint': 3, 'disp': True}, bounds=len(x0) * [[-.04, .04]])
        # 'L-BFGS-B'
        # bounds=len(x0)
        #                * [[-.04, .04]],
        Du, Dl = format_input(
            sol.x, self.bu.g, self.bu.g_p, self.bl.g, self.bl.g_p)
        self.bu.g.D = Du
        if sum(self.bl.g.dependent) != 0:
            self.bl.g.g_independent = self.bu.g
        self.bl.g.D = Dl
        # self.bu.g.internal_variables(self.bu.length, origin=self.bu.origin)
        # self.bl.g.internal_variables(self.bl.length, origin=self.bl.origin)
        self.calculate_x()

        self.bu.y = self.bu.g.x3(self.bu.x)
        self.bl.y = self.bl.g.x3(self.bl.x)
        print('sol', self.bu.g.D, self.bl.g.D)

    def formatted_residual(self, x0, format_input):
        [Au, Al] = format_input(x0, self.bu.g, self.bu.g_p, self.bl.g, self.bl.g_p)
        self.bu.l.concentrated_load = self.bu.l.external_load.copy()
        self.bl.l.concentrated_load = self.bl.l.external_load.copy()
        self.bu.l.concentrated_s = self.bu.l.external_s.copy()
        self.bl.l.concentrated_s = self.bl.l.external_s.copy()
        self.calculate_x()
        self.calculate_force()
        R = self.bu._residual(Au)
        if sum(self.bl.g.dependent) != 0:
            self.bl.g.g_independent = self.bu.g
        R += self.bl._residual(Al)
        if self.spars_s is not None:
            self.calculate_resultants()
        Ru = self.bu._residual(Au)
        if sum(self.bl.g.dependent) != 0:
            self.bl.g.g_independent = self.bu.g
        Rl = self.bl._residual(Al)
        R = .5*np.sqrt(Ru**2 + Rl**2)

        print('R', R, Rl, Ru)
        # BREAK
        return R

    def calculate_force(self):
        if self.coupled_forces:
            index = np.where(np.around(self.bu.s, decimals=7) ==
                             self.bu.l.concentrated_s[0])[0][0]
            x_u = self.bu.g.x1_grid[index]
            index = np.where(np.around(self.bl.s, decimals=7) ==
                             self.bl.l.concentrated_s[0])[0][0]
            x_l = self.bl.g.x1_grid[index]

            y_u = self.bu.g.x3(np.array([x_u]))[0]
            y_l = self.bl.g.x3(np.array([x_l]))[0]

            dx = x_u - x_l
            dy = y_u - y_l
            ds = np.sqrt(dx**2 + dy**2)

            Fx = self.bu.l.concentrated_magnitude[0]*dx/ds
            Fy = self.bu.l.concentrated_magnitude[0]*dy/ds
            self.bu.l.concentrated_load = [[-Fx, -Fy], ]
            self.bl.l.concentrated_load = [[Fx, Fy], ]
            self.bu.l.external_load = self.bu.l.concentrated_load.copy()
            self.bl.l.external_load = self.bl.l.concentrated_load.copy()

    def calculate_resultants(self):
        def M_f(s, *Rys):
            self.bl.l.concentrated_load = self.bl.l.external_load.copy()
            for Ry in Rys:
                self.bl.l.concentrated_load += [[0, Ry], ]
            self.bl.calculate_M()
            if self.bl.ignore_ends:
                return self.bl.M[1:-1]
            else:
                return self.bl.M
        self.bl.g.calculate_angles()

        M = (self.bl.p.young*self.bl.p.inertia)*(self.bl.g.rho - self.bl.g_p.rho)
        s = self.bl.s[1:-1]
        self.bu.l.concentrated_s = self.bu.l.external_s.copy() + self.spars_s
        self.bl.l.concentrated_s = self.bl.l.external_s.copy() + self.spars_s

        self.Rx = []
        self.Ry = []
        if self.spars_k is None:
            if self.bl.ignore_ends:
                M = M[1:-1]
                s = s[1:-1]
            popt, pcov = curve_fit(M_f, s, M, p0=[0]*len(self.spars_s))

            self.Rx += len(popt)*[0]  # -2.614  # sbeta*R
            self.Ry += list(popt)  # -3.634  # cbeta*R
        else:
            if len(self.spars_s) != 1:
                raise(NotImplementedError)
            else:
                i = 0
                # Current length for child
                index_u = np.where(np.around(self.bu.s, decimals=7) ==
                                   self.spars_s[0])[0][0]
                x_u = self.bu.g.x1_grid[index_u]
                index_l = np.where(np.around(self.bl.s, decimals=7) ==
                                   self.spars_s[0])[0][0]
                x_l = self.bl.g.x1_grid[index_l]

                y_u = self.bu.g.x3(np.array([x_u]))[0]
                y_l = self.bl.g.x3(np.array([x_l]))[0]

                dx = x_u - x_l
                dy = y_u - y_l
                ds = np.sqrt(dx**2 + dy**2)

                # Original length for parent
                x_u = self.bu.g_p.x1_grid[index_u]
                x_l = self.bl.g_p.x1_grid[index_l]

                y_u = self.bu.g_p.x3(np.array([x_u]))[0]
                y_l = self.bl.g_p.x3(np.array([x_l]))[0]

                dx = x_u - x_l
                dy = y_u - y_l
                ds0 = np.sqrt(dx**2 + dy**2)

                magnitude = self.spars_k[i]*(ds-ds0)

                self.Rx = [magnitude*dx/ds]
                self.Ry = [magnitude*dy/ds]

        self.bl.l.concentrated_load = self.bl.l.external_load.copy()
        self.bu.l.concentrated_load = self.bu.l.external_load.copy()
        for i in range(len(self.spars_s)):
            self.bl.l.concentrated_load += [[self.Rx[i], self.Ry[i]], ]
            self.bu.l.concentrated_load += [[-self.Rx[i], -self.Ry[i]], ]


def derivative(f, a, method='central', h=0.01):
    '''Compute the difference formula for f'(a) with step size h.

    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    method : string
        Difference formula: 'forward', 'backward' or 'central'
    h : number
        Step size in difference formula

    Returns
    -------
    float
        Difference formula:
            central: f(a+h) - f(a-h))/2h
            forward: f(a+h) - f(a))/h
            backward: f(a) - f(a-h))/h
    '''
    if method == 'central':
        return (f(a + h) - f(a - h))/(2*h)
    elif method == 'forward':
        return (f(a + h) - f(a))/h
    elif method == 'backward':
        return (f(a) - f(a - h))/h
    else:
        raise ValueError("Method must be 'central', 'forward' or 'backward'.")
