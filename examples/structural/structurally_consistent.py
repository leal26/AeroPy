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


def constraint_f(input):
    Au, Al = format_input(input, gu=a.bu.g, gu_p=a.bu.g_p, gl=a.bl.g, gl_p=a.bl.g_p)
    a.bu.g.D = Au
    a.bl.g.D = Al

    index_u = np.where(a.bu.s == a.spars_s[0])[0][0]
    index_l = np.where(a.bl.s == a.spars_s[0])[0][0]
    a.bu.g.calculate_x1(a.bu.s)
    a.bl.g.calculate_x1(a.bl.s)
    x_u = a.bu.g.x1_grid[index_u]
    x_l = a.bl.g.x1_grid[index_l]
    # x_u = a.bu.g.x1_grid[index]
    # x_l = a.bl.g.x1_grid[index]
    y_u = a.bu.g.x3(np.array([x_u]))[0]
    y_l = a.bl.g.x3(np.array([x_l]))[0]

    norm = math.sqrt((x_u-x_l)**2+(y_u-y_l)**2)
    s1 = (x_u - x_l)/norm
    s2 = (y_u - y_l)/norm
    xp_u = np.array([a.bu.g_p.x1_grid[index_u]])
    xp_l = np.array([a.bl.g_p.x1_grid[index_l]])
    delta = a.bu.g_p.x3(xp_u)[0] - a.bl.g_p.x3(xp_l)[0]
    a.bl.g.spar_directions = [[s1, s2]]
    print('Constraint', norm - delta, norm, delta)
    return norm - delta


def format_u(input, g=None, g_p=None):
    return list(input)


def format_input(input, gu=None, gu_p=None, gl=None, gl_p=None):
    _, _, n_u = g_upper._check_input([])

    Au = format_u(input[:n_u], gu, gu_p)
    Al = format_u(input[n_u:], gl, gl_p)
    return Au, Al


warnings.filterwarnings("ignore", category=RuntimeWarning)


psi_spars = [0.2]
m = len(psi_spars)

g_upper = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0.],
                                chord=[psi_spars[0], 1-psi_spars[0]],
                                color=['b', 'r'], N1=[1., 1.], N2=[1., 1.],
                                offset=.05, continuity='C2', free_end=True,
                                root_fixed=True)
g_lower = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0., 0., 0.],
                                chord=[psi_spars[0], 0.7, 0.1],
                                color=['b', 'r', 'g'], N1=[1., 1., 1.], N2=[1., 1., 1.],
                                offset=-.05, continuity='C2', free_end=True,
                                root_fixed=False, dependent=[True, False, False])

g_upper.calculate_s(N=[11, 9])
g_lower.calculate_s(N=[11, 8, 6])

# g_upper.D = [0.01, 0.02, 0.03, 0.04]
g_upper.D = [0, 0, 0, 0]
g_lower.g_independent = g_upper
# g_lower.D = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
g_lower.D = [0, 0, 0, 0, 0, 0]
g_upper.calculate_x1(g_upper.s)
g_lower.calculate_x1(g_lower.s)
print('D0', g_upper.cst[0].D, g_lower.cst[0].D)
print('D1', g_upper.cst[1].D, g_lower.cst[1].D)
print('D2', g_lower.cst[2].D)

# print('chords', a.bl.g.chord, a.bu.g.chord)
index_u = np.where(g_upper.s == psi_spars[0])[0][0]
index_l = np.where(g_lower.s == psi_spars[0])[0][0]
plt.figure()
plt.plot(g_upper.x1_grid, g_upper.x3(g_upper.x1_grid), 'b',
         label='Upper', lw=3)
plt.plot(g_lower.x1_grid, g_lower.x3(g_lower.x1_grid), 'b',
         label='Upper', lw=3)

plt.plot([g_upper.x1_grid[index_u], g_lower.x1_grid[index_l]], [g_upper.x3(
    np.array([g_upper.x1_grid[index_u]])), g_lower.x3(np.array([g_lower.x1_grid[index_l]]))], 'b', lw=3)
plt.show()
