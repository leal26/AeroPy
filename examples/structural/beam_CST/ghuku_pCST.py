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

# Ghuku paper
distributed_load = None
chord_parent = 0.4385
width = 0.0385
thickness = 0.00625
clamped_load = 20.9634
experiment = {0: [0., 0., 0.5492, 0.2000],
              138.321+clamped_load: [0., 0., 0.3746, 0.1753],
              211.896+clamped_load: [0., 0., 0.2475, 0.2864], }
# 287.433: [0., 0., 0.1287, 0.3582],
# 460.089: [0., 0., -0.1891, 0.6087],}
p = properties(young=200e9, dimensions=[width, thickness])
exp_F = [0.000, 20.9634, 83.7474, 159.2844, 232.8594, 308.3964, 382.9542, 459.4704, 481.0524]
exp_delta = [0, 0.0047, 0.0135, 0.03098, 0.0451, 0.0622, 0.0793, 0.0977, 0.103]
g0 = CoordinateSystem.polynomial(D=experiment[list(experiment.keys())[
    0]], chord=chord_parent, color='b')

colors = ['0.0', '0.3', '0.5', '0.7']
load_keys = list(experiment.keys())

g = CoordinateSystem.pCST(D=[-0.10248271, -0.10555923, -0.10863574, -0.19294965, -0.19987181], chord=[.4*chord_parent, .6*chord_parent], color=['b', 'r'],
                          N1=[1, 1], N2=[1, 1], continuity='C2', free_end=True, root_fixed=True)
g.calculate_s([21, 21])
s = g.s
for i in range(2, 3):  # len(load_keys)
    print('Load: ', load_keys[i])
    load = load_keys[i]
    l = loads(concentrated_load=[[0, -load]], load_s=[g.total_length],
              distributed_load=distributed_load)
    b = beam_chen(g, p, l, s, ignore_ends=False)
    # b.iterative_solver()
    # b._residual(b.g.D[:-1])
    b.parameterized_solver(format_input=format_input,
                           x0=b.g.D[:-1],
                           solver='gradient')

    if i == 0:
        plt.plot(b.x, b.y, colors[i], label='Parent', linestyle='-',
                 lw=3, zorder=1)
        gg = CoordinateSystem.polynomial(D=experiment[load_keys[i]], chord=chord_parent)
        gg.calculate_x1(s)
        gg.plot(label='Experiment: %.3f N' %
                load_keys[i], color=colors[i], scatter=True, marker="D")
    if i > 0:
        plt.plot(b.x, b.y, colors[i], label='Child: %.3f N' % load, linestyle='-',
                 lw=3, zorder=1)
        gg = CoordinateSystem.polynomial(D=experiment[load_keys[i]], chord=chord_parent)
        gg.calculate_x1(s)
        gg.plot(label='Experiment: %.3f N' %
                load_keys[i], color=colors[i], scatter=True, marker="D")
plt.legend()
plt.show()

plt.figure()
plt.plot(b.g.x1_grid[1:], b.M[1:], 'b', label='From forces')
M = (b.p.young*b.p.inertia)*(b.g.rho - b.g_p.rho)
plt.plot(b.g.x1_grid[1:], M[1:], 'r', label='From CST')
plt.legend()
plt.show()
