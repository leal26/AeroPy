import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def format_input(input):
    return [0, 0] + list(input)


abaqus_x = [0, 0.10398348, 0.20769666, 0.31089693, 0.4133383, 0.51480412,
            0.61510199, 0.71406609, 0.81157351, 0.90749156, 1.0017918]
abaqus_y = [0, 0.0025992838, 0.010444776, 0.023406046, 0.041411132,
            0.064284593, 0.091835514, 0.12384841, 0.16006343, 0.20029253, 0.24419737]

g = CoordinateSystem.polynomial(D=[0, 0, 0.25, 0], chord=1, color='b')
g_p = CoordinateSystem.polynomial(D=[0, 0, 0.25, 0], chord=1, color='k')

x = np.linspace(0, 1, 20)
s = np.zeros(len(x))
for i in range(len(x)):
    s[i] = g.arclength(x[i])[0]

p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[s[-1]])
# s = np.linspace(0, 1, 10)
b = beam_chen(g, p, l, s)
b.parameterized_solver(format_input=format_input, x0=b.g.D[2:])

g_p.calculate_x1(s)
g_p.plot(label='Parent')
plt.plot(b.x, b.y, '.5', label='Child: %.3f N' % -l.concentrated_load[0][-1], lw=3)
plt.scatter(abaqus_x, abaqus_y, c='.5', label='FEA: %.3f N' % -
            l.concentrated_load[0][-1], edgecolors='k', zorder=10, marker="^")

plt.legend()
plt.show()
