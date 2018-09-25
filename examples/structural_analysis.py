from aeropy.geometry.parametric import poly

import matplotlib.pyplot as plt
import numpy as np


curve = poly()


z1 = np.linspace(0, 1, 10)
z2 = curve.z2(z1)
t = curve.t(z1)
n = curve.n(z1)

g1 = curve.g1(z1)
g2 = curve.g2(z1)

h = .1
a = z1[4]
fprime = (curve.z2(a+h)-curve.z2(a-h))/h/2
tan = curve.z2(a) + fprime*(z1-a)

plt.figure()
plt.plot(z1, z2, label='geometry')
plt.quiver(z1, z2, g1[0], g1[1], angles='xy', label='CCS')
plt.quiver(z1, z2, g2[0], g2[1], angles='xy')
plt.legend()
plt.axis('equal')
plt.grid()
plt.show()
