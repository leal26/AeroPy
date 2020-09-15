import matplotlib.pyplot as plt
from scipy.optimize import fixed_point
import numpy as np
import pickle
import os
import math
import numpy as np
from numpy.linalg import inv
import warnings

from aeropy.structural.beam import beam_chen, coupled_beams
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem
from aeropy.geometry.airfoil import CST
from aeropy.CST_2D import calculate_c_baseline, calculate_psi_goal, calculate_spar_direction, S

warnings.filterwarnings("ignore")
g = CoordinateSystem.pCST(D=[.2, .3, .4, 1., 2.5, -.4], chord=[.2, .8], color=['b', 'r'],
                          N1=[1., 1], N2=[1, 1], continuity='C2')
print('lengths', g.cst[0].length, g.cst[1].length)


# g.x1_grid  = np.linspace(0, g.total_chord, 51)
g.calculate_s([11, 11])
# s = np.linspace(0, g.total_length, 21)
print('indexes', g.indexes, g.cst[0].indexes, g.cst[1].indexes)
g.calculate_x1(g.s)
print('D before', g.cst[0].D, g.cst[1].D)
g.D = [.2, .3, .4, 1., 2.5, -.4]
g.calculate_x1(g.s)
print('D after', g.cst[0].D, g.cst[1].D)
print('x1_grid', g.x1_grid)
dd = g.x3(g.x1_grid, diff='x11')

# dd2 = (1/g.cst[1].chord)*(-2*(g.nn+1)*g.cst[1].D[0] + 2*g.cst[1].D[1]*g.nn)
# dd1 = (1/g.cst[0].chord)*(2*g.nn*g.cst[0].D[-3] - 2*(g.cst[0].N1+g.nn)*g.cst[0].D[-2])
# print('DD', dd1, dd2, g.n, g.nn)
print('dd', dd)
print('chord', g.cst[0].chord, g.cst[1].chord, g.cst[0].chord + g.cst[1].chord)
print('offset', g.cst[0].offset_x, g.cst[1].offset_x)
plt.plot(g.x1_grid, dd)
plt.figure()
g.plot(label=['1', '2'])
plt.show()
