from aeropy.geometry.parametric import poly
import numpy as np
import math


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
