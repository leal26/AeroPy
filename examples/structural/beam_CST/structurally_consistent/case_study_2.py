import matplotlib.pyplot as plt
import warnings
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def constraint_f(input):
    A = format_input(input)
    g.D = A

    index = np.where(b.s == b.spars_s[0])[0][0]
    g.calculate_x1(b.s)
    x = g.x1_grid[index]
    y = g.x3(np.array([x]))

    return y


def format_input(input, g=None, g_p=None):
    # COnsidering BC for zero derivative at the root
    return list(input)
    # return input


warnings.filterwarnings("ignore", category=RuntimeWarning)

abaqus = np.loadtxt('case_study_2.csv', delimiter=',')

constraints = ({'type': 'eq', 'fun': constraint_f})
n = 2
p = 3
i = (n+1)*p
g = CoordinateSystem.pCST(D=i*[0], chord=[.2, .4, .4], color=['b', 'g', 'r'],
                          N1=[1, 1, 1], N2=[1, 1, 1], continuity='C1', free_end=True,
                          root_fixed=True)

p = properties()
l = loads(concentrated_load=[[0, -40]], load_s=[1], torque=[5], torque_s=[0.6])
g.calculate_s(N=[5, 9, 9])
b = beam_chen(g, p, l, s=None, calculate_resultant=True)
b.spars_s = [0.2]
b.parameterized_solver(format_input, x0=[0, 0, 0, 0, 0, 0, 0, 0], constraints=constraints)
# b.g.D = [2.3667520183461132e-07, 0.011413199268466722, 0.022841016335144255, 0.07336877500069414,
#          0.06464943617957472, 0.05578540569596733, 0.03693876321850288, 0.02798090460688548]
# b.g.calculate_x1(b.g.s)
# b.x = b.g.x1_grid
# b.y = b.g.x3(b.g.x1_grid)
# b.g.radius_curvature(b.g.x1_grid)
# b.calculate_resultants()

# Checking stuff
# print('check1', +b.g.cst[0].D[2] + b.g.cst[1].D[0] + b.g.cst[1].zetaT)
# print('check2', 2*b.g.cst[1].D[1] - (1+2)*b.g.cst[1].D[2], 2/3*b.g.cst[1].D[1])
# print('chord', b.g.cst[0].chord, b.g.cst[1].chord)
# print('offset', b.g.cst[0].offset_x, b.g.cst[1].offset_x)
# print('zetaT', b.g.cst[0].zetaT, b.g.cst[1].zetaT, b.g.zetaT)
# print('zetaL', b.g.cst[0].zetaL, b.g.cst[1].zetaL, b.g.zetaL)
# print('A0', b.g.A0)
print('loads', b.l.concentrated_load)
# print('D', b.g.cst[0].D, b.g.cst[1].D)
plt.figure()
plt.plot(b.g.x1_grid[1:], b.M[1:], 'b', label='From forces')
print('M1', b.M[1:])

M = (b.p.young*b.p.inertia)*(b.g.rho - b.g_p.rho)
print('M2', M[1:])
for i in range(len(M)):
    print(i, b.g.x1_grid[i], M[i], b.M[i])
print('indexes', b.g.cst[0].indexes, b.g.cst[1].indexes, b.g.cst[2].indexes)
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
plt.plot(b.g.x1_grid, b.g.x3(b.g.x1_grid, diff='x11'), 'r',
         label='dd', lw=3)
plt.legend()
plt.show()
