import math
import numpy as np
from scipy import integrate


class poly():
    """class for a polynomial function
    """

    def __init__(self, a=[-math.sqrt(3)/3, math.sqrt(3)/3]):
        self.a = a

    def z2(self, z1):
        """ z2 (checked)"""
        return(self.a[0]*z1**3 + self.a[1]*z1**2)

    def z2z1(self, z1):
        """ dz2 / dz1 (checked)"""
        return(3*self.a[0]*z1**2 + 2*self.a[1]*z1)

    def z2z11(self, z1):
        """ d^2z2 / dz1^2 (checked)"""
        return(6*self.a[0]*z1 + 2*self.a[1])

    def x1z1(self, z1):
        """ dx1/ dz1 (checked)"""
        return(np.sqrt(1+(self.z2z1(z1))**2))

    def z1x1(self, z1):
        """ dx1 / dz1 (checked)"""
        return(1.0/self.x1z1(z1))

    def z2x1(self, z1):
        """ dx2 / dx1 (checked)"""
        return(self.z2z1(z1)*self.z1x1(z1))

    def z2x11(self, z1):
        """ d^2z2 / dx1^2 """
        return(self.z2z11(z1)*(self.z1x1(z1))**2 +
               self.z2z1(z1)*self.z1x11(z1))

    def x1z11(self, z1):
        """ d^2x1 / dz1^2 """
        return(self.z1x1(z1)*self.z2z1(z1)*self.z2z11(z1))

    def z1x11(self, z1):
        """ d^2z1 / dz1^2 """
        return(self.z1x1(z1)**3*self.z2z1(z1)*self.z2z11(z1))

    def x1(self, z_final):
        """ d^2z2 / dz1^2 (checked)"""
        output, err = integrate.quad(self.x1z1, 0, z_final)
        return(output)

    def g1(self, z1, normalize=True):
        """Tangent vector r (checked)"""
        try:
            output = np.array([1, self.z2z1(z1)])
        except(ValueError):
            output = np.array([np.ones(len(z1)), self.z2z1(z1)])
        if normalize:
            output = output*self.z1x1(z1)
        return(output)

    def g2(self, z1, diff=False, normalize=True):
        """Normal vector (checked)"""
        if diff:
            try:
                output = np.array([- self.z2z11(z1), np.zeros(len(z1))])
            except(TypeError):
                output = np.array([- self.z2z11(z1), 0])
        else:
            try:
                output = np.array([- self.z2z1(z1), np.ones(len(z1))])
            except(TypeError):
                output = np.array([- self.z2z1(z1), 1])
        if normalize:
            return(output*self.z1x1(z1))
        else:
            return(output)

    def neutral_line(self, z1):
        """ Position along neutral line"""
        return(np.array([z1, self.a[0]*z1**3 + self.a[1]*z1**2]))

    def r(self, z1, x2):
        """ Position anywhere along shell considering shell thickness """
        return(self.neutral_line(z1) + x2*self.g2(z1))

    def rx1(self, z1):
        """ Position anywhere along shell considering shell thickness """
        output = np.array([self.z1x1(z1), self.z2x1(z1)])
        return(output)

    def rx11(self, z1):
        """ Position anywhere along shell considering shell thickness """
        return(np.array([self.z1x11(z1), self.z2x11(z1)]))

    def g11(self, z1, x2=0):
        g1 = self.g1(z1)
        return(np.inner(g1, g1))

    def g22(self, z1):
        g2 = self.g2(z1)
        return(np.inner(g2, g2))


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
