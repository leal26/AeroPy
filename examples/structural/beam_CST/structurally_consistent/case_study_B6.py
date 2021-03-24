import matplotlib.pyplot as plt
import warnings
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem

from scipy.interpolate import interp1d

def rmse(x1, y1, x2, y2):
    print('x1', x1)
    print('x2', x2)
    kind = "cubic"
    x1[x1<0] = 0
    x2[x2<0] = 0
    x1 = np.delete(x1, 21)
    y1 = np.delete(y1, 21)
    if max(x1) > max(x2):
        f = interp1d(x1, y1, kind=kind)
        predictions = y2
        targets = f(x2)
    else:
        f = interp1d(x2, y2, kind=kind)
        predictions = y1
        targets = f(x1)
    return np.sqrt(np.mean((predictions-targets)**2))

def format_input(input, g=None, g_p=None):
    # COnsidering BC for zero derivative at the root
    return list(input)
    # return input


warnings.filterwarnings("ignore", category=RuntimeWarning)

abaqus = np.loadtxt('case_study_B6.csv', delimiter=',')

D_n4 = [0.288197, 0.37918722, 0.21903394, 0.25925831, 0.24259769, 0.27278782,
        0.06528661, 0.07283469, 0.06692653, 0.0559131, -0.01268226, -0.00134043,
        -0.08743061, -0.0013288]

D_n3 = [3.06132404e-01, 3.31838944e-01, 2.44074856e-01, 2.30082260e-01,
        2.71512524e-01, 4.98187521e-02, 8.59903075e-02, 4.89862447e-02,
        1.22835553e-04, -7.05141464e-02, -1.95302773e-02]

D_n2 = [0.3208729, 0.28457915, 0.2431048, 0.27140013, 0.06470603, 0.06222951,
        -0.02810757, -0.04412038]

# Db_n3 = [0.30645304, 0.33092568, 0.2454664, 0.23047635, 0.27141484, 0.04881688,
#          0.09722178, 0.02404006]

# g = CoordinateSystem.pCST(D=Db_n3,
#                           chord=[.2, .8], color=['b', 'r', 'g'],
#                           N1=[.5, 1], N2=[1, 1], continuity='C2',
#                           free_end=True, rigid_LE=True)

g = CoordinateSystem.pCST(D=D_n2,
                          chord=[.2, .7, .1], color=['b', 'r', 'g'],
                          N1=[.5, 1, 1], N2=[1, 1, 1], continuity='C2',
                          free_end=True, rigid_LE=True, root_fixed=False)

p = properties()
l = loads(concentrated_load=[[0, -10]], load_s=[g.cst[0].length + g.cst[1].length])
g.calculate_s(N=[21, 21, 21])
b = beam_chen(g, p, l, s=None, ignore_ends=False)
b.g.calculate_x1(b.g.s)
b.x = b.g.x1_grid
b.y = b.g.x3(b.x)
b.g.radius_curvature(b.g.x1_grid)
b.calculate_M()
_, _, n_u = g._check_input([])
print('D1', b.g.cst[0].D)
print('D2', b.g.cst[1].D)
print('D3', b.g.cst[2].D)
# b.parameterized_solver(format_input, x0=np.array([0.11055139, 0.05173331, 0.09983107]))
# b.parameterized_solver(format_input, x0=np.array(
#     [0.10688066710180918, 0.04881687846018838, 0.0972217791492235]))
b.parameterized_solver(format_input, x0=np.array(
    [0.08817500333333335, 0.06470603, 0.06222951, -0.02810757]), solver='lm')
# b.parameterized_solver(format_input, x0=np.zeros(n_u))
# b._residual([0.08817500333333335, 0.06470603, 0.06222951, -0.02810757])
# 0.06470603, 0.06222951, -0.02810757
print('Residual', b.R)
print('r', b.r, len(b.r))
print('s', b.g.s, len(b.g.s))
print('D1', b.g.cst[0].D)
print('D2', b.g.cst[1].D)
print('D3', b.g.cst[2].D)
# print('chord', b.g.cst[0].chord, b.g.cst[1].chord, b.g.cst[2].chord)
# print('offset_x', b.g.cst[0].offset_x, b.g.cst[1].offset_x, b.g.cst[2].offset_x)
# print('zetaL', b.g.cst[0].zetaL, b.g.cst[1].zetaL, b.g.cst[2].zetaL)
# print('zetaT', b.g.cst[0].zetaT, b.g.cst[1].zetaT, b.g.cst[2].zetaT)
# print('zL', b.g.cst[0].zetaL*b.g.cst[0].chord, b.g.cst[1].zetaL *
#       b.g.cst[1].chord, b.g.cst[2].zetaL*b.g.cst[2].chord)
# print('zT', b.g.cst[0].zetaT*b.g.cst[0].chord, b.g.cst[1].zetaT *
#       b.g.cst[1].chord, b.g.cst[2].zetaT*b.g.cst[2].chord)
# print('D1', b.g.cst[0].D)
# print('D2', b.g.cst[1].D)
# print('D3', b.g.cst[2].D)
# print('loads', b.l.concentrated_load)

import matplotlib
matplotlib.rcParams.update({'font.size': 14})

print(rmse(b.g.x1_grid, b.g.x3(b.g.x1_grid), abaqus[0, :], abaqus[1, :]))

plt.figure()
plt.plot(b.s, b.M, '.5', lw=3, label='From forces', clip_on=False)

M = (b.p.young*b.p.inertia)*(b.g.rho - b.g_p.rho)
plt.scatter(b.s, M[21:], c='.5', edgecolors='k', zorder=20, marker='D', label='From CST', clip_on=False)
plt.xlim([0, max(b.s)])
plt.ylabel("Units (N.m)")
plt.xlabel("s (m)")
plt.legend()

plt.figure()
plt.plot(b.g_p.x1_grid, b.g_p.x3(b.g_p.x1_grid), 'k',
         label='Upper Parent', lw=3, clip_on=False)
plt.plot(b.g.x1_grid, b.g.x3(b.g.x1_grid), c='.5',
         label='Upper Child', lw=3, clip_on=False)
plt.scatter(abaqus[0, :], abaqus[1, :], c='.5',
            label='FEA', edgecolors='k', zorder=10, marker="^", clip_on=False)

plt.figure()
plt.plot(b.g.x1_grid, b.g.x3(b.g.x1_grid, diff='x1'), 'b',
         label='d', lw=3)
plt.plot(b.g_p.x1_grid, b.g_p.x3(b.g_p.x1_grid, diff='x1'), 'b--',
         label='d', lw=3)
plt.plot(b.g.x1_grid, b.g.x3(b.g.x1_grid, diff='x11'), 'r',
         label='dd', lw=3)
plt.plot(b.g_p.x1_grid, b.g_p.x3(b.g_p.x1_grid, diff='x11'), 'r--',
         label='dd', lw=3)
plt.legend()
plt.show()
