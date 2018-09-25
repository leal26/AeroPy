from aeropy.geometry.parametric import poly

import matplotlib.pyplot as plt
import numpy as np


curve = poly()


z1 = np.linspace(0, 1, 10)
z2 = curve.z2(z1)
g1 = curve.g1(z1)
g2 = curve.g2(z1)
rx1 = curve.rx1(z1)
rx11 = curve.rx11(z1)

plt.figure()
plt.plot(z1, z2, label='geometry')
plt.quiver(z1, z2, g1[0], g1[1], angles='xy', label='CCS')
plt.quiver(z1, z2, g2[0], g2[1], angles='xy')
plt.quiver(z1, z2, rx1[0], rx1[1], angles='xy', color='g')
plt.quiver(z1, z2, rx11[0], rx11[1], angles='xy', color='g')
plt.legend()
plt.axis('equal')
plt.grid()
plt.show()
