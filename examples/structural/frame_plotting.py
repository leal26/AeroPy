from aeropy.geometry.parametric import poly
from aeropy.geometry.frame import frame

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import numpy as np


curve = poly()
z1 = np.linspace(0, 1)
z2 = curve.z2(z1)

colors = ['r', 'g', 'k', 'b', 'm']
plt.figure()
plt.plot(z1, z2, label='geometry')
r1 = curve.r(z1, diff=[0, 1, 0, 0])
r2 = curve.r(z1, diff=[0, 0, 1, 0])
plt.plot(r1[0], r1[1], label='a1')
plt.plot(r2[0], r2[1], label='a2')
z1 = np.linspace(.1, 1, 6)
i = 0
x1 = curve.x1(z1)

for x2 in np.linspace(-.3, .3, 5):
    g1 = curve.g(1, z1, x2)
    g2 = curve.g(2, z1, x2)

    r = curve.r(z1, x2)
    Q = plt.quiver(r[0], r[1], g1[0], g1[1], angles='xy', color=colors[i],
                   label='$x^2$=%.2f' % x2, scale_units='xy')
    plt.quiver(r[0], r[1], g2[0], g2[1], angles='xy', color=colors[i],
               scale_units='xy')
    i += 1
plt.legend(loc='best')
plt.axis('equal')
plt.xlabel('$z_1$')
plt.ylabel('$z_2$')
plt.show()
