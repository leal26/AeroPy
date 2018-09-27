from aeropy.geometry.parametric import poly, frame

import matplotlib.pyplot as plt
import numpy as np


curve = poly()
z1 = np.linspace(0, 1)
z2 = curve.z2(z1)

ccs = frame(curve=curve, frame='Bishop')
[T, N, B] = ccs.frenet_serret()
g1 = curve.g1(ccs.z1)
g2 = curve.g2(ccs.z1)

plt.figure()
plt.plot(z1, z2, label='geometry')
plt.quiver(ccs.z1, ccs.z2, g1[0], g1[1], angles='xy', label='Truth')
plt.quiver(ccs.z1, ccs.z2, g2[0], g2[1], angles='xy')
plt.quiver(ccs.z1, ccs.z2, T[0], T[1], angles='xy', color='g',
           label='Bishop')
plt.quiver(ccs.z1, ccs.z2, N[0], N[1], angles='xy', color='g')
plt.legend()
plt.axis('equal')
plt.show()
