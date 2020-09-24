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
# g = CoordinateSystem.pCST(D=[.2, .3, .3, .4, .4, 1., -.4], chord=[.2, .8], color=['b', 'r'],
#                           N1=[1., 1], N2=[1, 1], continuity='C1', free_end=True)
g = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0., 0., 0.], chord=[.2, .7, .1], color=['b', 'r', 'g'],
                          N1=[1., 1, 1], N2=[1, 1, 1], continuity='C2', free_end=True, root_fixed=True)
g.calculate_s([11, 11, 5])
# g.D = [0, 0, 0, 0, 0, 0]
g.D = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
g.calculate_x1(g.s)
print('chord', g.cst[0].chord, g.cst[1].chord, g.cst[2].chord)
print('offset_x', g.cst[0].offset_x, g.cst[1].offset_x, g.cst[2].offset_x)
print('zetaL', g.cst[0].zetaL, g.cst[1].zetaL, g.cst[2].zetaL)
print('zetaT', g.cst[0].zetaT, g.cst[1].zetaT, g.cst[2].zetaT)
print('zL', g.cst[0].zetaL*g.cst[0].chord, g.cst[1].zetaL *
      g.cst[1].chord, g.cst[2].zetaL*g.cst[2].chord)
print('zT', g.cst[0].zetaT*g.cst[0].chord, g.cst[1].zetaT *
      g.cst[1].chord, g.cst[2].zetaT*g.cst[2].chord)
print('D1', g.cst[0].D)
print('D2', g.cst[1].D)
print('D3', g.cst[2].D)
dd = g.x3(g.x1_grid, diff='x11')
plt.plot(g.x1_grid, dd)
plt.figure()
g.plot(label=['1', '2', '3'])
plt.show()
