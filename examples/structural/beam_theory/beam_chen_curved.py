import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem

abaqus_x = [0, 0.050000664, 0.100002654, 0.150012404, 0.200025707, 0.250053465,
            0.300086498, 0.350139558, 0.400199026, 0.450283021, 0.50037396,
            0.550492764, 0.600618601, 0.650774479, 0.700936854, 0.751130462,
            0.801329434, 0.851559758, 0.901793659, 0.952057898, 1.002323747]
abaqus_y = [0, 0.001222659, 0.002391184, 0.006008377, 0.009577097, 0.015600263,
            0.021580685, 0.030021306, 0.038424935, 0.049294505, 0.060132839,
            0.073442832, 0.086727314, 0.102489136, 0.118231162, 0.136456192,
            0.154667079, 0.175366625, 0.196057677, 0.219242975, 0.242425412]

g = CoordinateSystem.polynomial(D=[0, 0, 0.25, 0], chord = 1, color = 'b')

p = properties()
l = loads(concentrated_load = [[0, -1]], load_s = [1])
s = np.linspace(0, 1, 10)
x0 = s
b = beam_chen(g, p, l, s)
b.iterative_solver()

x = b.x
y = g.x3(x0) + b.y

plt.plot(x, y, 'b', label='Chen', linestyle = '--', lw = 3)
plt.scatter(abaqus_x, abaqus_y, c='g', label='FEA', edgecolors='k', zorder = 10)
plt.legend()
plt.show()
