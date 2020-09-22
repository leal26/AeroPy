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

g.D = [-0.00020088858978972732, 0.000200894085615972, -0.00020088905395165244, -
       0.00020088854026064262, -0.0002008882375582824, -0.0002008879895743138]
g.calculate_x1(g.s)
print('chord', g.cst[0].chord, g.cst[1].chord, g.cst[2].chord)
print('offset_x', g.cst[0].offset_x, g.cst[1].offset_x, g.cst[2].offset_x)
print('D1', g.cst[0].D)
print('D2', g.cst[1].D)
print('D3', g.cst[2].D)
dd = g.x3(g.x1_grid, diff='x11')
plt.plot(g.x1_grid, dd)
plt.figure()
g.plot(label=['1', '2', '3'])
plt.show()
