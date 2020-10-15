import matplotlib.pyplot as plt
import warnings
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def format_input(input, g=None, g_p=None):
    # COnsidering BC for zero derivative at the root
    return list(input)
    # return input


warnings.filterwarnings("ignore", category=RuntimeWarning)

abaqus = np.loadtxt('case_study_C1.csv', delimiter=',')

D_n4 = [0.288197, 0.37918722, 0.21903394, 0.25925831, 0.24259769, 0.27278782,
        0.06528661, 0.07283469, 0.06692653, 0.0559131, -0.01268226, -0.00134043,
        -0.08743061, -0.0013288]

D_n3 = [3.06132404e-01, 3.31838944e-01, 2.44074856e-01, 2.30082260e-01,
        2.71512524e-01, 4.98187521e-02, 8.59903075e-02, 4.89862447e-02,
        1.22835553e-04, -7.05141464e-02, -1.95302773e-02]

D_n2 = [0.3208729, 0.28457915, 0.2431048, 0.27140013, 0.06470603, 0.06222951,
        -0.02810757, -0.04412038]
g = CoordinateSystem.pCST(D=D_n4,
                          chord=[.2, .7, .1], color=['b', 'r', 'g'],
                          N1=[.5, 1, 1], N2=[1, 1, 1], continuity='C2',
                          free_end=True, rigid_LE=True)

p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[g.cst[0].length + g.cst[1].length])
g.calculate_s(N=[21, 11, 11])
b = beam_chen(g, p, l, s=None, ignore_ends=False)
# b.g.calculate_x1(b.g.s)
# b.x = b.g.x1_grid
# b.y = b.g.x3(b.x)
# b.g.radius_curvature(b.g.x1_grid)
# b.calculate_M()
_, _, n_u = g._check_input([])
b.parameterized_solver(format_input, x0=np.zeros(n_u))
print('D1', b.g.cst[0].D)
print('D2', b.g.cst[1].D)
print('D3', b.g.cst[2].D)
print('loads', b.l.concentrated_load)

plt.figure()
print('x', np.shape(b.g.x1_grid[21:]))
print('M', np.shape(b.M[:]))
print('rho', np.shape(b.g.rho))
plt.plot(b.g.x1_grid[21:], b.M[:], 'b', label='From forces')

M = (b.p.young*b.p.inertia)*(b.g.rho - b.g_p.rho)
plt.plot(b.g.x1_grid[1:], M[1:], 'r', label='From CST')
plt.legend()

plt.figure()
plt.plot(b.g_p.x1_grid, b.g_p.x3(b.g_p.x1_grid), 'b',
         label='Upper Parent', lw=3)
plt.plot(b.g.x1_grid, b.g.x3(b.g.x1_grid), c='.5',
         label='Upper Child', lw=3)
plt.scatter(abaqus[0, :], abaqus[1, :], c='.5',
            label='FEA', edgecolors='k', zorder=10, marker="^")

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
