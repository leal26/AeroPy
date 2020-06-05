import aeropy
import math

import numpy as np
from scipy.integrate import quad
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

    def update_chord(self, length_target = None, bounds = None):
        def f(c_c):
            length_current, err = self.g_c.arclength(c_c)
            return abs(length_target - length_current)
        if length_target is None:
            length_target= self.arc_length
        if bounds is None:
            self.g_c.chord = minimize(f, self.g_c.chord).x[0]
        else:
            self.g_c.chord = minimize(f, self.g_c.chord,
                                               method='L-BFGS-B',
                                               bounds = bounds).x[0]
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
            self.W = np.trapz(np.multiply(self.load,u), self.g_c.x1_grid)

    def residual(self):
        self.R = self.U - self.W

    def minimum_potential(self, x0=[0,0]):
        def to_optimize(x):
            self.g_c.D = [0,0] + list(x)
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
        bounds = np.array(((-0.01,0.01),)*len(x0))

        res = minimize(to_optimize, x0)
        self.g_c.D = [0,0] + list(res.x)
        self.work()
        self.free_energy()
        self.strain_energy()
        self.residual()
        return(res.x, res.fun)

class beam_chen():
    def __init__(self, geometry, properties, load, s):
        self.g = geometry
        self.p = properties
        self.l = load
        self.s = s
        self.length = self.g.arclength()

    def G(self, x):
        c = 1/self.p.young/self.p.inertia
        # print(x)
        # return c*(quad(self.M, 0, x)[0])
        c = self.g.chord
        return self.l.concentrated_load[0][-1]/self.p.young/self.p.inertia*(c*x-x**2/2)

    def M(self, x):
        c = self.g.chord
        for i in range(len(self.l.concentrated_s)):
            concentrated_x_i = self._x(self.l.concentrated_s[i])
            return self.l.concentrated_load[0][-1]*(concentrated_x_i*-x)

    def s_to_x(self):
        self.x = np.zeros(len(self.s))
        for i in range(len(self.s)):
            self.x[i] = self._x(self.s[i])

    def _x(self, s):
        def _to_minimize(l):
            def _to_integrate(x):
                den = np.sqrt(1-self.G(x)**2)
                if np.isnan(den):
                    return 100
                else:
                    return 1/den
            l = l[0]
            current_L = quad(_to_integrate, 0, l)[0]
            return abs(s-current_L)
        return minimize(_to_minimize, s, method = 'Nelder-Mead',).x[0]

    def find_deflection(self):
        def _to_integrate(x):
            G = self.G(x)
            den = np.sqrt(1-G**2)
            if np.isnan(den):
                return 100
            else:
                return G/den

        self.y = []
        for x_i in self.x:
            y_i = current_L = quad(_to_integrate, 0, x_i)[0]
            self.y.append(y_i)
