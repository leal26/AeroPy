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
n = 2
p = 3
i = n*p+2
g = CoordinateSystem.pCST(D=np.linspace(0.01, 0.01*i, i), chord=[.2, .7, .1], color=['b', 'r', 'g'],
                          N1=[.5, 1, 1], N2=[1, 1, 1], continuity='C2', free_end=True, rigid_LE=True)
# g = CoordinateSystem.pCST(D=[3.06132404e-01, 3.31838944e-01, 2.44074856e-01,
#                              2.30082260e-01, 4.98187521e-02, 8.59903075e-02,
#                              4.89862447e-02, 1.22835553e-04, -7.05141464e-02,
#                              -1.95302773e-02,  2.71512524e-01], chord=[.2, .7, .1], color=['b', 'r', 'g'],
#                           N1=[.5, 1, 1], N2=[1, 1, 1], continuity='C2', free_end=True, root_fixed=True,
#                           rigid_LE=True)

g.calculate_s([21, 11, 11])
g.D = [0.07416666666666666, 0.05, 0.06, 0.07]
# g.D = [0.01, 0.02, 0.03]
g.calculate_x1(g.s)
print('x', g.x1_grid)

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
dd = g.x3(g.x1_grid, diff='x1')
plt.plot(g.x1_grid, dd)
plt.figure()
g.plot(label=['1', '2', '3'])
plt.show()
