import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def format_input(input, g=None, g_p=None):
    # COnsidering BC for zero derivative at the root
    return list(input) + [-input[0]]


abaqus = np.loadtxt('output_A7.csv', delimiter=',')

g = CoordinateSystem.CST(D=[0.20000000000000004, 0.2500000000000002, 0.2583333333333329,
                            0.22500000000000048, 0.1500000000000001, -0.20000000000000004], chord=1, color='b', N1=1, N2=1)
g.calculate_s(21)

g.calculate_x1(g.s)
p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[g.s[-1]])

b = beam_chen(g, p, l, s=None, ignore_ends=True)
# b.parameterized_solver(format_input, x0=b.g.D[:-1])
x0 = list(np.array(b.g.D[:-1]) + np.array([0.00571429, 0.005,
                                           0.00428571, 0.00357143, 0.00285714]))
b.parameterized_solver(format_input, x0=x0)
# b._residual(x0)

# data = np.array([b.x, b.y]).T
# np.savetxt('input_A7.csv', data, delimiter=',')

plt.figure()
plt.plot(b.g_p.x1_grid, b.g_p.x3(b.g_p.x1_grid), 'b',
         label='Upper Parent', lw=3)
plt.plot(b.x, b.y, '.5', label='Child: %.3f N' % -l.concentrated_load[0][1],
         lw=3, zorder=0)
plt.scatter(abaqus[0, :], abaqus[1, :], c='.5',
            label='FEA', edgecolors='k', zorder=10, marker="^")

plt.legend()

plt.figure()
plt.plot(b.g.x1_grid, b.g.rho)

plt.figure()
plt.plot(b.g.x1_grid[1:], b.M[1:], 'b', label='From forces')
M = (b.p.young*b.p.inertia)*(b.g.rho - b.g_p.rho)
plt.plot(b.g.x1_grid[1:], M[1:], 'r', label='From CST')
plt.legend()

plt.figure()
plt.plot(b.g.x1_grid[1:], b.g.x3(b.g.x1_grid, diff='x11')[1:], 'b', label='x11')
plt.plot(b.g.x1_grid[1:], b.g.rho[1:], '-r', label='rho')
plt.legend()
plt.show()
