from scipy.optimize import approx_fprime
import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem

from aeropy.CST_2D import dxi_u, dxi_l, ddxi_u, ddxi_l
from aeropy.geometry.airfoil import CST

def format_input(input):
    return input[2:4]

g = CoordinateSystem.CST(D=[0, 0, 0.25, 0], chord = 1, color = 'b')

x = np.linspace(0.0, 1, 21)
# s = np.zeros(len(x))
# for i in range(len(x)):
#     s[i] = g.arclength(x[i])[0]

g.x1_grid = x
print(g.D[:-2])
d_original = -dxi_l(x, g.D[:-1], 0)
dd_original = -ddxi_l(x, g.D[:-1])
d_new = dxi_u(x, g.D[:-1], 0)
dd_new = ddxi_u(x, g.D[:-1], 0)

y = CST(x, 1, deltasz=0, Au=g.D[:-1])
print(y)
print('old: ', d_original)
print('new: ', d_new)

plt.figure()
plt.plot(x, y, 'k', label='CST')
plt.plot(x, d_original, label='Original')
plt.plot(x, dd_original, 'r--', label='DD')
plt.scatter(x, d_new, label='New')
plt.scatter(x, dd_new, label='DDNew')
plt.legend()
plt.show()
