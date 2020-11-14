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

n = 2
p = 2
i = n*p+1
g_upper = CoordinateSystem.pCST(D=i*[0.],
                                chord=[psi_spars[0], 1-psi_spars[0]],
                                color=['b', 'r'], N1=[1., 1.], N2=[1., 1.],
                                offset=.05, continuity='C2', free_end=True,
                                root_fixed=True)
n = 2
p = 3
i = n*p+1
g_lower = CoordinateSystem.pCST(D=i*[0.],
                                chord=[psi_spars[0], 0.7, 0.1],
                                color=['b', 'r', 'g'], N1=[1., 1., 1.], N2=[1., 1., 1.],
                                offset=-.05, continuity='C2', free_end=True,
                                root_fixed=True)

g_upper.calculate_s(N=[11, 9])
g_lower.calculate_s(N=[11, 8, 6])
p_upper = properties()
p_lower = properties()
l_upper = loads(concentrated_load=[[-np.sqrt(2)/2, -np.sqrt(2)/2]], load_s=[1])
l_lower = loads(concentrated_load=[[np.sqrt(2)/2, np.sqrt(2)/2]], load_s=[1-0.1])

a = coupled_beams(g_upper, g_lower, p_upper, p_lower, l_upper, l_lower, None,
                  None)

a.calculate_x()

# a.formatted_residual(format_input=format_input, x0=[
#                      0.00200144, 0.00350643, 0.00255035, 0.00226923] + [-0.00219846, - 0.00313221, - 0.00193564, - 0.00191324])
# a.formatted_residual(format_input=format_input, x0=[
#                      0.00200144, 0.00350643, 0.00255035, 0.00226923, 0.00183999] + [-0.00219846, - 0.00313221, - 0.00193564, - 0.00191324, - 0.00127513])
# a.formatted_residual(format_input=format_input, x0=list(
# g_upper.D[:-1]) + list(g_lower.D[:1]) + list(g_lower.D[2:-1]))
_, _, n_u = g_upper._check_input([])
_, _, n_l = g_lower._check_input([])
a.parameterized_solver(format_input=format_input, x0=list(
    g_upper.D[:n_u]) + list(g_lower.D[:n_l]), solver='lm')
# input = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
# [Au, Al] = format_input(input, a.bu.g, a.bu.g_p, a.bl.g, a.bl.g_p)
# a.calculate_x()
# a.calculate_force()
# R = a.bu._residual(Au) + a.bl._residual(Al)
# if a.spars_s is not None:
#     a.calculate_resultants()
# R = a.bu._residual(Au) + a.bl._residual(Al)
print('upper', a.bu.g.D)
print('upper 1', a.bu.g.cst[0].D)
print('upper 2', a.bu.g.cst[1].D)
print('lower', a.bl.g.D)
print('lower 1', a.bl.g.cst[0].D)
print('lower 2', a.bl.g.cst[1].D)
print('lower 3', a.bl.g.cst[2].D)
print('loads', a.bl.l.concentrated_load, a.bu.l.concentrated_load)
plt.figure()
plt.plot(a.bu.g.x1_grid[1:], a.bu.M[1:], 'b', label='Upper')
plt.plot(a.bl.g.x1_grid[1:], a.bl.M[1:], 'r', label='Lower')

Ml = (a.bl.p.young*a.bl.p.inertia)*(a.bl.g.rho - a.bl.g_p.rho)
Mu = (a.bu.p.young*a.bu.p.inertia)*(a.bu.g.rho - a.bu.g_p.rho)
plt.plot(a.bu.g.x1_grid[1:], Mu[1:], '--b', label='Upper')
plt.plot(a.bl.g.x1_grid[1:], Ml[1:], '--r', label='Lower')
plt.legend()

# print('chords', a.bl.g.chord, a.bu.g.chord)
plt.figure()
plt.plot(a.bu.g_p.x1_grid, a.bu.g_p.x3(a.bu.g_p.x1_grid), 'b',
         label='Upper Parent', lw=3)
plt.plot(a.bu.g.x1_grid, a.bu.g.x3(a.bu.g.x1_grid), c='.5',
         label='Upper Child: %.3f N' % -l_upper.concentrated_load[0][-1], lw=3)
plt.plot(a.bl.g_p.x1_grid, a.bl.g_p.x3(a.bl.g_p.x1_grid), 'b', linestyle='dashed',
         label='Lower Parent', lw=3)
plt.plot(a.bl.g.x1_grid, a.bl.g.x3(a.bl.g.x1_grid), '.5', linestyle='dashed',
         label='Lower Child: %.3f N' % -l_upper.concentrated_load[0][-1], lw=3)

upper = np.loadtxt('case_study_4_upper.csv', delimiter=',')
lower = np.loadtxt('case_study_4_lower.csv', delimiter=',')
plt.scatter(upper[0, :], upper[1, :], c='.5', label='Abaqus', edgecolors='k',
            zorder=10, marker="^")
plt.scatter(lower[0, :], lower[1, :], c='.5', edgecolors='k', zorder=10,
            marker="^")

# x = [a.bu.g.chord*a.bl.g.spar_psi_upper[0], a.bl.g.chord*a.bl.g.spar_psi[0]]
# y = [a.bu.g.chord*a.bl.g.spar_xi_upper[0], a.bl.g.chord*a.bl.g.spar_xi[0]]
# dx = x[1]-x[0]
# dy = y[1]-y[0]
# norm = math.sqrt(dx**2+dy**2)
# print('spar direction', a.bl.g.spar_directions)
# print('actual direction', dx/norm, dy/norm)
# plt.plot(x, y, c='g', label='spars', lw=3)
# plt.arrow(x[0], y[0], -a.bl.g.spar_directions[0][0]*a.bl.g.delta_P[0],
#           -a.bl.g.spar_directions[0][1]*a.bl.g.delta_P[0])
# print(a.bl.g.delta_P[0])
plt.legend()
# plt.gca().set_aspect('equal', adjustable='box')
plt.show()
