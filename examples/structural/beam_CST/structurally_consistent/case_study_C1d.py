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

abaqus = np.loadtxt('case_study_C1d.csv', delimiter=',')

n = 2
p = 2
i = n*p+2
D1 = [0]*(n+1)
m = i - n - 2
D2 = list(np.linspace(0.2, 0.2*m, m))
x0 = D1 + D2[:-1]
g = CoordinateSystem.pCST(D=D1+D2, chord=[.2, .8], color=['b', 'r'],
                          N1=[1, 1], N2=[1, 1], continuity='C2', free_end=True, root_fixed=True)

p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[g.cst[0].length + g.cst[1].length])
g.calculate_s(N=[11, 31])
b = beam_chen(g, p, l, s=None, ignore_ends=False)
b.g.calculate_x1(b.g.s)
b.x = b.g.x1_grid
b.y = b.g.x3(b.x)
b.g.radius_curvature(b.g.x1_grid)
b.calculate_M()
_, _, n_u = g._check_input([])
print('D1', b.g.cst[0].D)
print('D2', b.g.cst[1].D)
# print('D3', b.g.cst[2].D)
# b.parameterized_solver(format_input, x0=np.array([0.11055139, 0.05173331, 0.09983107]))
# b.parameterized_solver(format_input, x0=np.array(
#     [0.10688066710180918, 04r3et5r.04881687846018838, 0.0972217791492235]))
# b.parameterized_solver(format_input, x0=np.array(
#     [0.08817500333333335, 0.06470603, 0.06222951, -0.02810757]))
# b.parameterized_solver(format_input, x0=x0)
# x0 = [0.0016609031927165981, 0.0014583095210324206, 0.0014721447209376166, -
#          0.0015928952913609152, 0.10242922551427547, 0.20275245669550482, -0.20183072908959696]
# x0 = [0.0016609031927165981, 0.0014583095210324206, 0.0014721447209376166,
#          0.10242922551427547, 0.20275245669550482, -0.20183072908959696]
# x0 = [0.0160733147016186, -0.002225084058028139, -0.0020682982405123436, 0.18571910365412023]
# b._residual(x0)
plt.figure()
# try:
b.parameterized_solver(format_input, x0=x0)
# except:
#     plt.legend()
#     plt.ylim([-1, 1])
#     plt.show()
print('Residual', b.R)
print('r', b.r, len(b.r))
print('s', b.g.s, len(b.g.s))
print('D1', b.g.cst[0].D)
print('D2', b.g.cst[1].D)
print('lengths', b.g.cst[0].length, b.g.cst[1].length)
# print('D3', b.g.cst[2].D)
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

plt.figure()
print('x', np.shape(b.g.x1_grid[:]))
print('M', np.shape(b.M[:]))
print('rho', np.shape(b.g.rho))
plt.plot(b.g.x1_grid[:], b.M[:], 'b', label='From forces')

M = (b.p.young*b.p.inertia)*(b.g.rho - b.g_p.rho)
plt.plot(b.g.x1_grid, M, 'r', label='From CST')
plt.plot(b.g.x1_grid[:], (b.p.young*b.p.inertia)*b.r, 'g', label='residual')
plt.legend()

plt.figure()
plt.plot(b.g_p.x1_grid, b.g_p.x3(b.g_p.x1_grid), 'b',
         label='Upper Parent', lw=3)
plt.plot(b.g.x1_grid, b.g.x3(b.g.x1_grid), c='.5',
         label='Upper Child', lw=3)
plt.scatter(abaqus[:, 0], abaqus[:, 1], c='.5',
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
