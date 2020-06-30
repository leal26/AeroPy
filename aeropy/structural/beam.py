import aeropy
import math
import copy

import numpy as np
from scipy.integrate import quad, trapz
from scipy.optimize import minimize


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
        bounds = np.array(((-0.01, 0.01),)*len(x0))

        res = minimize(to_optimize, x0)
        self.g_c.D = [0, 0] + list(res.x)
        self.work()
        self.free_energy()
        self.strain_energy()
        self.residual()
        return(res.x, res.fun)


class beam_chen():
    def __init__(self, geometry, properties, load, s):
        self.g = copy.deepcopy(geometry)
        self.g_p = copy.deepcopy(geometry)
        self.g_p.calculate_x1(s)
        # self.Rho = self.g_p.x3(self.g_p.x1_grid, diff='theta11')
        self.g_p.radius_curvature(self.g_p.x1_grid)
        self.p = properties
        self.l = load
        self.s = s
        self.length = self.g.arclength()[0]

    def calculate_M(self):
        self.M = np.zeros(len(self.x))
        for i in range(len(self.M)):
            self.M[i] = self._M(self.x[i], self.s[i], self.y[i])

    def _M(self, x, s, y=None):
        M_i = 0
        if self.l.concentrated_load is not None:
            if self.l.follower:
                raise(NotImplementedError)
            else:
                for i in range(len(self.l.concentrated_s)):
                    index = np.where(self.s == self.l.concentrated_s[i])[0][0]
                    M_i += self.l.concentrated_load[i][0]*(self.y[index]-y)
                    M_i += self.l.concentrated_load[i][1]*(self.x[index]-x)

        if self.l.distributed_load is not None:
            index = np.where(self.s == s)[0][0]

            if not self.l.follower:
                # M_x + M_y
                M_i -= trapz(self.l.distributed_load(self.s[index:])
                             * (self.x[index:]-x), self.s[index:])
            else:
                w = self.l.distributed_load(self.s[index:])
                M_x = w*self.cos[index:]*(self.x[index:]-x)
                print(self.sin)
                M_y = w*self.sin[index:]*(self.y[index:]-y)
                M_i -= trapz(M_x + M_y, self.s[index:])
                # print('cos', self.cos[index:])
                # print('sin', self.sin[index:])
                # print('x', self.x[index:])
                # print('y', self.y[index:])
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

    def calculate_angles(self):
        self.cos = self.g.x1(self.g.x1_grid, 'theta1')
        self.sin = self.g.x3(self.g.x1_grid, 'theta1')

    def calculate_deflection(self):
        self.y = np.zeros(len(self.x))
        for i in range(len(self.x)):
            dydx = self.G[:i+1]/(1-self.G[:i+1]**2)
            y_i = trapz(dydx, self.x[:i+1])
            self.y[i] = y_i

    # def calculate_residual(self):
    #     self.r = np.zeros(len(self.x))
    #     for i in range(len(self.x)):
    #         rhs = self.G[i]/(1-self.G[i]**2)
    #         lhs = self.g.x3(self.x[i], diff='x1')
    #         self.r[i] = lhs - rhs
    #     self.R = np.linalg.norm(self.r)

    def calculate_residual(self, ignore_ends=False):
        self.r = np.zeros(len(self.x))
        for i in range(len(self.x)):
            # rhs = self.g.x3(self.x[i], diff='theta11') - self.Rho[i]
            rhs = self.g.rho[i] - self.g_p.rho[i]
            lhs = self.M[i]/self.p.young/self.p.inertia
            self.r[i] = np.abs(lhs - rhs)
        # self.R = np.linalg.norm(self.r)
        if ignore_ends:
            self.R = trapz(self.r[1:-1], self.s[1:-1])
        else:
            self.R = trapz(self.r, self.s)
        # print('r', self.r)
        if np.isnan(self.R):
            self.R = 100
        # print('r', self.r)
        print('R: ', self.R)
    # def calculate_Pho(self):
    #     self.g.calculate_x1(self.theta1, bounds = self.g_p.bounds)
    #     self.g.basis()
    #     self.g.basis(diff = 'theta')
    #     self.g.metric_tensor()
    #     self.g.metric_tensor(diff = 'theta')
    #     self.g.curvature_tensor()
    #     self.Pho = self.B[0][0]

    def iterative_solver(self):
        # Calculated undeformed properties
        # self.calculate_Pho()
        # Calculate Moment for undeformed
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

    def parameterized_solver(self, format_input=None, x0=None,
                             ignore_ends=False):
        def formatted_residual(A):
            A = format_input(A)
            return self._residual(A, ignore_ends=ignore_ends)

        sol = minimize(formatted_residual, x0, method='SLSQP')
        self.g.D = format_input(sol.x)
        self.y = self.g.x3(self.x)
        print('sol', self.g.D, sol.fun)
        return

    def _residual(self, A, ignore_ends=False):
        self.g.D = A
        self.g.calculate_x1(self.s)
        self.g.chord = self.g.x1_grid[-1]
        self.x = self.g.x1_grid
        self.y = self.g.x3(self.x)
        if self.l.follower:
            self.calculate_angles()
        self.calculate_M()
        # self.calculate_G()
        # self.calculate_x()
        # self.calculate_M()
        self.g.radius_curvature(self.g.x1_grid)
        self.calculate_residual(ignore_ends)
        return self.R
