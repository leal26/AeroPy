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
    a.bu.g.internal_variables(arc_upper)
    a.bl.g.internal_variables(arc_upper)

    index = np.where(a.bl.s == a.spars_s[0])[0][0]

    x_u = a.bu.g.x1_grid[index]
    x_l = a.bl.g.x1_grid[index]
    y_u = a.bu.g.x3(np.array([x_u]))[0]
    y_l = a.bl.g.x3(np.array([x_l]))[0]

    norm = math.sqrt((x_u-x_l)**2+(y_u-y_l)**2)
    s1 = (x_u - x_l)/norm
    s2 = (y_u - y_l)/norm
    x_p = np.array([a.bu.g_p.x1_grid[index]])
    delta = a.bu.g_p.x3(x_p)[0] - a.bl.g_p.x3(x_p)[0]
    a.bl.g.spar_directions = [[s1, s2]]
    return norm - delta


def constraint_ff(Au, Al):
    a.bu.g.D = Au
    a.bl.g.D = Al
    a.bu.g.internal_variables(arc_upper)
    a.bl.g.internal_variables(arc_upper)

    index = np.where(a.bl.s == a.spars_s[0])[0][0]

    x_u = a.bu.g.x1_grid[index]
    x_l = a.bl.g.x1_grid[index]
    y_u = a.bu.g.x3(np.array([x_u]))[0]
    y_l = a.bl.g.x3(np.array([x_l]))[0]

    norm = math.sqrt((x_u-x_l)**2+(y_u-y_l)**2)
    s1 = (x_u - x_l)/norm
    s2 = (y_u - y_l)/norm
    x_p = np.array([a.bu.g_p.x1_grid[index]])
    delta = a.bu.g_p.x3(x_p)[0] - a.bl.g_p.x3(x_p)[0]
    a.bl.g.spar_directions = [[s1, s2]]
    return norm - delta


def format_u(input, g=None, g_p=None):
    # COnsidering BC for zero derivative at the root
    def free_end(An_in):
        g.D = list(input) + [An_in, -input[0]]
        g.internal_variables(arc_upper)
        den = (1+(-An_in + g.zetaT)**2)**(1.5)
        An_out = (2*n*Cn1 - den*(g.chord/g_p.chord)*dd_p)/(2*N1+2*n)
        return An_out

    # A5
    n = len(g_p.D) - 2
    Pn = g_p.D[-2]
    Pn1 = g_p.D[-3]
    Cn1 = input[-1]
    N1 = g.N1
    dd_p = (2*n*Pn1-2*(N1+n)*Pn)/(1+(-Pn + g_p.zetaT)**2)**(1.5)
    An = float(fixed_point(free_end, g_p.D[-2]))
    return list(input) + [An, -input[0]]
    # return input


def format_input(input, gu=None, gu_p=None, gl=None, gl_p=None):
    Au = format_u(input[:int(len(input)/2)], gu, gu_p)
    Al = format_u(input[int(len(input)/2):], gl, gl_p)
    # print('Constraint', constraint_ff(Au, Al))
    return Au, Al


warnings.filterwarnings("ignore", category=RuntimeWarning)

upper = np.loadtxt('coupled_beam_upper.csv', delimiter=',')
lower = np.loadtxt('coupled_beam_lower.csv', delimiter=',')

psi_spars = [0.2]
m = len(psi_spars)

g_upper = CoordinateSystem.CST(D=[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., ], chord=1, color='b', N1=1, N2=1,
                               offset=.05)
g_lower = CoordinateSystem.CST(D=[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., ], chord=1, color='b', N1=1, N2=1,
                               offset=-.05)
g_p = CoordinateSystem.CST(D=[0., 0., 0., 0., 0., 0., 0., 0., 0.,
                              0., 0.], chord=1, color='k', N1=1, N2=1)

s_upper = np.linspace(0, g_p.arclength(np.array([1])), 21)
s_lower = np.linspace(0, g_p.arclength(np.array([1])), 21)
p_upper = properties()
p_lower = properties()
l_upper = loads(concentrated_load=[[-np.sqrt(2)/2, -np.sqrt(2)/2]], load_s=[1])
l_lower = loads(concentrated_load=[[np.sqrt(2)/2, np.sqrt(2)/2]], load_s=[1-0.1])
# l_upper = loads(concentrated_load=[[0, -1]], load_s=[1])
# l_lower = loads(concentrated_load=[[0, 1]], load_s=[1])
arc_upper = 1.0
arc_lower = 1.0
arc_spar = g_lower.arclength(np.array([psi_spars[0]]))

a = coupled_beams(g_upper, g_lower, p_upper, p_lower, l_upper, l_lower, s_upper,
                  s_lower, ignore_ends=True, spars_s=psi_spars)
a.calculate_x()
constraints = ({'type': 'eq', 'fun': constraint_f})
# a.formatted_residual(format_input=format_input, x0=[
#                      0.00200144, 0.00350643, 0.00255035, 0.00226923] + [-0.00219846, - 0.00313221, - 0.00193564, - 0.00191324])
# a.formatted_residual(format_input=format_input, x0=[
#                      0.00200144, 0.00350643, 0.00255035, 0.00226923, 0.00183999] + [-0.00219846, - 0.00313221, - 0.00193564, - 0.00191324, - 0.00127513])
# a.formatted_residual(format_input=format_input, x0=list(
# g_upper.D[:-1]) + list(g_lower.D[:1]) + list(g_lower.D[2:-1]))
a.parameterized_solver(format_input=format_input, x0=list(
    g_upper.D[:-2]) + list(g_lower.D[:-2]), constraints=constraints)
print(a.bu.g.D, a.bl.g.D)
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
index = np.where(a.bl.s == a.spars_s[0])[0][0]
plt.figure()
plt.plot(a.bu.g_p.x1_grid, a.bu.g_p.x3(a.bu.g_p.x1_grid), 'b',
         label='Upper Parent', lw=3)
plt.plot(a.bu.g.x1_grid, a.bu.g.x3(a.bu.g.x1_grid), c='.5',
         label='Upper Child: %.3f N' % -l_upper.concentrated_load[0][-1], lw=3)
plt.plot(a.bl.g_p.x1_grid, a.bl.g_p.x3(a.bl.g_p.x1_grid), 'b', linestyle='dashed',
         label='Lower Parent', lw=3)
plt.plot(a.bl.g.x1_grid, a.bl.g.x3(a.bl.g.x1_grid), '.5', linestyle='dashed',
         label='Lower Child: %.3f N' % -l_upper.concentrated_load[0][-1], lw=3)
plt.plot([a.bu.g_p.x1_grid[index], a.bl.g_p.x1_grid[index]], [a.bu.g_p.x3(
    a.bu.g_p.x1_grid[index]), a.bl.g_p.x3(a.bl.g_p.x1_grid[index])], 'b', lw=3)
plt.plot([a.bu.g.x1_grid[index], a.bl.g.x1_grid[index]], [a.bu.g.x3(
    a.bu.g.x1_grid[index]), a.bl.g.x3(a.bl.g.x1_grid[index])], '.5', lw=3)
plt.scatter(upper[0, :], upper[1, :], c='.5', label='FEA: %.3f N' % -
            l_upper.concentrated_load[0][-1], edgecolors='k', zorder=10, marker="^")
plt.scatter(lower[0, :], lower[1, :], c='.5', edgecolors='k', zorder=10, marker="^")
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
