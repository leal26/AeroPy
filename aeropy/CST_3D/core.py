import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from math import factorial


def piecewise_linear(x, y):
    """Returns piecewise linear function based x-y pairs provided"""
    return interp1d(x, y)


class BernsteinPolynomial:
    """Implements Bernstein polynomial of arbitrary order"""

    def __init__(self, order, coefficients=None):
        self._order = order
        if coefficients is not None:
            self._coefficients = coefficients
        else:
            self._coefficients = n_coeff*[1.]
        self._K = self._calculate_k(order+1)

    @staticmethod
    def _calculate_k(number):
        K = number*[0]
        n = number-1
        for r in range(n+1):
            K[r] = factorial(n)/(factorial(r)*factorial(n-r))

        return K

    def set_coefficients(self, coeff):
        self._coefficients = coeff

    def __call__(self, param1, param2=None):

        A = self._coefficients
        K = self._K
        n = self._order
        try:
            F = np.zeros(param1.shape)
        except AttributeError:
            F = 0.

        for r in range(n+1):
            F += (try_as_func(A[r], param2)*K[r]*np.power(param1, r) *
                  np.power(1.-param1, n-r))

        return F


def try_as_func(f, x, y=None):
    """
    if f is function of x and y, returns f(x, y)
    if f is function, returns f(x)
    if f is a float, returns f
    """
    if y is not None:
        try:
            fx = f(x, y)
        except TypeError:
            try:
                fx = f(x)
            except TypeError:
                fx = f
    else:
        try:
            fx = f(x)
        except TypeError:
            fx = f

    return fx


class CST3D:
    """Implements general Class/Shape Transformation function

    :param location: changes part location in global
    :param rotation: changes part rotation in global
    :param ref: changes reference in parameterized domain
    :param sx and sy:
    :param nx and ny:
    :param XYZ:
    :param xshear and yshear:"""

    def __init__(self, **params):
        self.location = params.get('location', (0., 0., 0.))
        self.rotation = params.get('rotation', (0., 0., 0.))
        self.ref = params.get('ref', (0., 0., 0.))
        self.sx = params.get('sx', 1.)
        self.nx = params.get('nx', (0., 0.))
        self.sy = params.get('sy', 1.)
        self.ny = params.get('ny', (0.5, 1.))
        self.XYZ = params.get('XYZ', (1., 1., 1.))
        self.xshear = params.get('xshear', 0.)
        self.zshear = params.get('zshear', 0.)
        self.twist = params.get('twist', 0.)
        self.TE_thickness = params.get('TE_thickness', 0.)

    def __call__(self, psi, eta):
        # calculate x, y, and z from the 3D CST equation
        x_l = self._cst_x(psi, eta)
        y_l = self._cst_y(eta)
        z_l = self._cst_z(psi, eta)

        # rotate
        twist = try_as_func(self.twist, eta)*np.pi/180.
        x_lr = np.cos(twist)*x_l+np.sin(twist)*z_l
        z_lr = -np.sin(twist)*x_l+np.cos(twist)*z_l

        # shear
        x_lrs = x_lr + try_as_func(self.xshear, eta)
        z_lrs = z_lr + try_as_func(self.zshear, eta)

        # rotate from local CST coordinate system to global
        x_g, y_g, z_g = self._local_to_global(x_lrs, y_l, z_lrs)

        return x_g, y_g, z_g
        # return x_l, y_l, z_l

    def inverse(self, x_g, y_g, z_g):
        x_lrs, y_l, z_lrs = self._global_to_local(x_g, y_g, z_g)
        import matplotlib.pyplot as plt

        eta = np.array(self._inverse_y(y_l))

        # Correction for in case slightly above 1 or below 0
        eta[eta < 0] = 0
        eta[eta > 1] = 1
        # shear
        x_lr = x_lrs - try_as_func(self.xshear, eta)
        z_lr = z_lrs - try_as_func(self.zshear, eta)

        # rotate
        twist = try_as_func(self.twist, eta)*np.pi/180.
        x_l = np.cos(twist)*x_lr-np.sin(twist)*z_lr
        z_l = np.sin(twist)*x_lr+np.cos(twist)*z_lr

        psi = self._inverse_x(x_l, eta)
        # plt.figure()
        # plt.scatter(psi, eta)
        # plt.title('check')
        # plt.show()
        return psi, eta

    def _cst_x(self, psi, eta):
        psi0 = self.ref[0]
        X = self.XYZ[0]

        sc_y = try_as_func(self.sy, eta)*self._cy(eta)
        x = sc_y*(psi-psi0)*X  # +try_as_func(self.xshear, eta)

        return x

    def _cst_y(self, eta):
        eta0 = self.ref[1]
        Y = self.XYZ[1]

        y = (eta-eta0)*Y

        return y

    def _cst_z(self, psi, eta):
        zeta0 = self.ref[2]
        Z = self.XYZ[2]

        sc_x = try_as_func(self.sx, psi, eta)*self._cx(psi, eta)
        sc_y = try_as_func(self.sy, eta)*self._cy(eta)
        z = (sc_x*sc_y-zeta0 + psi*self.TE_thickness(eta))*Z  # +try_as_func(self.zshear, eta)

        return z

    def _cx(self, psi, eta):
        # class function in x
        nx1 = try_as_func(self.nx[0], eta)
        nx2 = try_as_func(self.nx[1], eta)
        # kx = self._k(nx1, nx2)
        c = np.power(psi, nx1)*np.power(1.-psi, nx2)

        return c

    def _cy(self, eta):
        # class function in y
        ny1, ny2 = self.ny
        ky = self._k(ny1, ny2)
        c = ky*np.power(eta, ny1)*np.power(1.-eta, ny2)

        return c

    def _k(self, N1, N2):
        # constant that normalizes class function based on
        # class coefficients
        k = np.power(N1+N2, N1+N2)/(np.power(N1, N1)*np.power(N2, N2))
        return k

    def _inverse_y(self, y_l):

        eta = y_l/self.XYZ[1]+self.ref[1]

        return eta

    def _inverse_x(self, x_l, eta):
        psi0, eta0, zeta0 = self.ref
        X, Y, Z = self.XYZ

        sc_y = try_as_func(self.sy, eta)*self._cy(eta)
        # psi = (x_l-try_as_func(self.xshear, eta))/(sc_y*X)+psi0
        psi = (x_l)/(sc_y*X)+psi0

        return psi

    def _local_to_global(self, x_l, y_l, z_l):
        q = Euler2Quat(self.rotation)
        q = NormalizeQuaternion(q)
        x_g, y_g, z_g = Body2Fixed((x_l, y_l, z_l), q)
        dx, dy, dz = self.location

        x_g += dx
        y_g += dy
        z_g += dz

        return x_g, y_g, z_g

    def _global_to_local(self, x_g, y_g, z_g):
        q = Euler2Quat(self.rotation)
        q = NormalizeQuaternion(q)

        dx, dy, dz = self.location

        x_l, y_l, z_l = Fixed2Body((x_g-dx, y_g-dy, z_g-dz), q)

        return x_l, y_l, z_l


def intersection(wing, fuselage, psi_w, eta_w0):
    x_intersect = np.zeros(len(psi_w))
    y_intersect = np.zeros(len(psi_w))
    z_intersect = np.zeros(len(psi_w))
    for i, psi in enumerate(psi_w):
        eta = eta_w0
        # x_w, y_w, z_w = 1., 1., 1.,
        # x_f, y_f, z_f = 0., 0., 0.,
        error = 1.  # np.sqrt((x_w-x_f)**2+(y_w-y_f)**2+(z_w-z_f)**2)
        count = 0
        while np.abs(error) > 1.e-13 and count < 30:
            count += 1
            # print("psi, eta", psi, eta)
            x_w, y_w, z_w = wing(psi, eta)
            # print("psi, eta", psi, eta)
            # print("wing point:", x_w, y_w, z_w)
            # x_history.append(x_w)
            # y_history.append(y_w)
            # z_history.append(z_w)
            psi_f, eta_f = fuselage.inverse(x_w, y_w, z_w)
            x_f, y_f, z_f = fuselage(psi_f, eta_f)
            # print("fuselage point:", x_f, y_f, z_f)
            # x_history.append(x_f)
            # y_history.append(y_f)
            # z_history.append(z_f)
            toss, eta = wing.inverse(x_f, y_f, z_f)
            # print(psi, "psi_error", psi-toss)
            error = np.sqrt((x_w-x_f)**2+(y_w-y_f)**2+(z_w-z_f)**2)
            # print(psi, "error", error)
            # print("error: ", error)
            # input("Press Enter to continue....")
        # history = np.array([x_history, y_history, z_history]).T
        # np.savetxt("history"+str(i)+".csv", history, delimiter=',')
        x_w, y_w, z_w = wing(psi, eta)

        x_intersect[i] = x_w
        y_intersect[i] = y_w
        z_intersect[i] = z_w

    return x_intersect, y_intersect, z_intersect


def Euler2Quat(euler_angles):
    phi, theta, psi = [a*np.pi/180. for a in euler_angles]
    C0 = np.cos(phi/2.0)
    C1 = np.cos(theta/2.0)
    C2 = np.cos(psi/2.0)
    S0 = np.sin(phi/2.0)
    S1 = np.sin(theta/2.0)
    S2 = np.sin(psi/2.0)
    return [C0*C1*C2 + S0*S1*S2,
            S0*C1*C2 - C0*S1*S2,
            C0*S1*C2 + S0*C1*S2,
            C0*C1*S2 - S0*S1*C2]


def NormalizeQuaternion(q):
    return [x*(1.5 - 0.5*(q[0]**2 + q[1]**2 + q[2]**2 + q[3]**2)) for x in q]


def local_to_assembly(x, y, z, q):
    q_mat = [[q[1]**2 - q[2]**2 - q[3]**2,
              2.0*(q[1]*q[2]),
              2.0*(q[1]*q[3])],
             [2.0*(q[1]*q[2]),
              q[2]**2 - q[1]**2 - q[3]**2,
              2.0*(q[2]*q[3])],
             [2.0*(q[1]*q[3]),
              2.0*(q[3]*q[2]),
              q[3]**2 - q[1]**2 - q[2]**2]]
    return [q_mat[0][0]*x + q_mat[0][1]*y + q_mat[0][2]*z,
            q_mat[1][0]*x + q_mat[1][1]*y + q_mat[1][2]*z,
            q_mat[2][0]*x + q_mat[2][1]*y + q_mat[2][2]*z]


def Body2Fixed(vector, e):
    """Transforms a vector in the body-fixed frame to the earth-fixed."""

    v = [0., vector[0], vector[1], vector[2]]
    v_f = QuatMult(e, QuatMult(v, [e[0], -e[1], -e[2], -e[3]]))
    # v_f = e*(v*e.inverse())

    return v_f[1], v_f[2], v_f[3]


def Fixed2Body(vector, e):
    """Transforms a vector in the earth-fixed frame to the body-fixed."""

    v = [0., vector[0], vector[1], vector[2]]
    v_b = QuatMult([e[0], -e[1], -e[2], -e[3]], QuatMult(v, e))
    # v_b = e.inverse()*(v*e)

    return v_b[1], v_b[2], v_b[3]


def QuatMult(q1, q2):
    """Performs a quaternion product."""

    q10, q1x, q1y, q1z = q1
    q20, q2x, q2y, q2z = q2

    Q0 = (q10*q20 -
          q1x*q2x -
          q1y*q2y -
          q1z*q2z)
    Qx = (q10*q2x +
          q1x*q20 +
          q1y*q2z -
          q1z*q2y)
    Qy = (q10*q2y -
          q1x*q2z +
          q1y*q20 +
          q1z*q2x)
    Qz = (q10*q2z +
          q1x*q2y -
          q1y*q2x +
          q1z*q20)

    return Q0, Qx, Qy, Qz
