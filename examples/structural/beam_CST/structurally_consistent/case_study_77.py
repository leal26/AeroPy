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
    a.bl.g.g_independent = a.bu.g
    a.bl.g.D = Al
    # length_diff = []
    # for i in range(a.bl.g.p):
    #     if a.bl.g.dependent[i]:
    #         current_length = a.bl.g.cst[i].arclength(a.bl.g.cst[i].chord)
    #         target_length = a.bl.g.cst[i].length
    #         length_diff.append(target_length-current_length)
    # offset_x = self.cst[j].offset_x + self.cst[j].chord
    a.bl.g.cst[0].chord = a.bl.g.spar_x[0]
    current_length = a.bl.g.cst[0].arclength(a.bl.g.cst[0].chord)
    target_length = a.bl.g.cst[0].length
    length_diff = target_length-current_length
    print('C', length_diff, target_length, current_length, a.bl.g.spar_x[0], a.bl.g.delta_P)
    return np.array([length_diff])


def format_u(input, g=None, g_p=None):
    # [0.00856]
    return list(input)


def format_l(input, g=None, g_p=None):
    return list(input)


def format_input(input, gu=None, gu_p=None, gl=None, gl_p=None):
    _, _, n_u = g_upper._check_input([])
    n_u = n_u  # - np.count_nonzero(gl.dependent)

    Au = format_u(input[:n_u], gu, gu_p)
    Al = format_l(input[n_u:], gl, gl_p)
    return Au, Al


warnings.filterwarnings("ignore", category=RuntimeWarning)


psi_spars = [0.2]
chords = []
for i in range(len(psi_spars)):
    if i == 0:
        chords.append(psi_spars[i])
    else:
        chords.append(psi_spars[i] - psi_spars[i-1])
chords.append(1-psi_spars[-1])

m = len(psi_spars)

n = 2
p = 2
i = n*p+1

g_upper = CoordinateSystem.pCST(D=[ii*0.0 for ii in range(i)],
                                chord=chords, color=['b', 'r', 'g', 'm'],
                                N1=len(chords)*[1.], N2=len(chords)*[1.],
                                offset=.05, continuity='C2', free_end=True,
                                root_fixed=True)
n = 2
p = 3
i = n*p+1
chords = [psi_spars[0], 0.7, 0.1]
g_lower = CoordinateSystem.pCST(D=[ii*0.0 for ii in range(i)],
                                chord=chords, color=['b', 'r', 'g', 'm'],
                                N1=len(chords)*[1.], N2=len(chords)*[1.],
                                offset=-.05, continuity='C2', free_end=True,
                                root_fixed=True,
                                dependent=[True, False, False])

g_upper.calculate_s(N=[11, 9])
g_lower.calculate_s(N=[11, 8, 101])
p_upper = properties()
p_lower = properties()
l_upper = loads(concentrated_load=[[-1*np.sqrt(2)/2, -1*np.sqrt(2)/2]], load_s=[1])
l_lower = loads(concentrated_load=[[1*np.sqrt(2)/2, 1*np.sqrt(2)/2]], load_s=[1-0.1])


a = coupled_beams(g_upper, g_lower, p_upper, p_lower, l_upper, l_lower, None,
                  None, ignore_ends=True, spars_s=psi_spars)

a.calculate_x()
constraints = ({'type': 'eq', 'fun': constraint_f})
# a.formatted_residual(format_input=format_input, x0=[
#                      0.00200144, 0.00350643, 0.00255035, 0.00226923] + [-0.00219846, - 0.00313221, - 0.00193564, - 0.00191324])
# a.formatted_residual(format_input=format_input, x0=[
#                      0.00200144, 0.00350643, 0.00255035, 0.00226923, 0.00183999] + [-0.00219846, - 0.00313221, - 0.00193564, - 0.00191324, - 0.00127513])
# a.formatted_residual(format_input=format_input, x0=list(
# g_upper.D[:-1]) + list(g_lower.D[:1]) + list(g_lower.D[2:-1]))
_, _, n_u = g_upper._check_input([])
_, _, n_l = g_lower._check_input([])

a.parameterized_solver(format_input=format_input, x0=np.zeros(
    n_u+n_l), solver='lm')

# x0 = [0.01, 0.02, 0.03, 0.04, 0.05,
#       1.38314286e-02, 1.33855583e-02, 1.53252454e-02,
#       6.09392225e-03, 2.28507716e-06]
# x0 = [6.25098946e-02, 6.39223659e-02, 1.16135430e-01, 1.30855583e-02, 3.15241173e-03,
#       6.09392225e-03, 6.01280457e-03, 1.42014104e-02, 3.90711975e-03, 2.28507716e-06]
# Current best fit
# x0 = [-2.74723157e-04, -6.78225692e-03, 9.57503164e-02, 7.85615379e-04, -
#       3.13934468e-03, -3.13079422e-02, -1.97834229e-02, -4.94540315e-05]
# x0 = [1.77792950e-05, -5.92255597e-01,  1.67623961e-02, -
#       1.57451120e-03, -3.32890461e-01,  2.34412685e-01, -8.66890589e-02]
# 1 N
# x0 = [-1.35639619e-04, -1.41483687e-04,  2.01141408e-03, -2.64840663e-04,
#       -4.59880256e-04, -1.47928883e-03, -9.86870379e-04, 4.63518281e-08]

# 100 N
# x0 = [1.77792950e-05, -5.92255597e-02,  1.67623961e-02, -
#       1.57451120e-03, -3.32890461e-01,  2.34412685e-01, -8.66890589e-02]
# Au, Al = format_input(x0, gl=a.bl.g)
# a.update(Au, Al)
# x0 = [-2.48380157e-04,  6.37431602e-04,  1.13267966e-05,  2.26737726e-03,
#       2.02524662e-03,  1.00096256e-04, -1.20017531e-03, -1.10681663e-03,
#       -1.76396286e-03, -1.63801966e-03, -1.27157595e-03,  1.50637021e-04,
#       3.57213301e-04]

# fitted
# x0 = [9.79649849e-05, -8.91443987e-04,  1.63805021e-03, -2.05489855e-04,
#       -4.29515355e-04, -1.37183124e-03, -9.98246075e-04,  9.20083949e-05]
#
# R = a.formatted_residual(x0, format_input)
# a.calculate_x()

# print('R', R)
print('upper', a.bu.g.D)
print('upper 1', a.bu.g.cst[0].D)
print('upper 2', a.bu.g.cst[1].D)
print('lower', a.bl.g.D)
print('lower 1', a.bl.g.cst[0].D)
print('lower 2', a.bl.g.cst[1].D)
print('lower 3', a.bl.g.cst[2].D)
print('zetaT', a.bu.g.cst[0].zetaT, a.bu.g.cst[1].zetaT, a.bu.g.cst[0].chord*a.bu.g.cst[0].zetaT,
      a.bu.g.cst[1].chord*a.bu.g.cst[1].zetaT, a.bu.g.cst[0].chord, a.bu.g.cst[1].chord)
print('loads', a.bl.l.concentrated_load, a.bu.l.concentrated_load)
xi = np.array([a.bu.g.cst[0].chord])
print('TEST', xi, a.bu.g.x3(xi))
# print('lengths', a.bl.g.cst[0].length, a.bl.g.cst[1].length,
#       a.bl.g.cst[2].length, a.bl.g.cst[3].length)

plt.figure()
plt.plot(a.bu.g.x1_grid[1:], a.bu.M[1:], 'b', label='Upper (F)')
plt.plot(a.bl.g.x1_grid[1:], a.bl.M[1:], 'r', label='Lower (F)')

Ml = (a.bl.p.young*a.bl.p.inertia)*(a.bl.g.rho - a.bl.g_p.rho)
Mu = (a.bu.p.young*a.bu.p.inertia)*(a.bu.g.rho - a.bu.g_p.rho)
plt.plot(a.bu.g.x1_grid[1:], Mu[1:], '--b', label='Upper (CST)')
plt.plot(a.bl.g.x1_grid[1:], Ml[1:], '--r', label='Lower (CST)')
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


for i in range(len(a.spars_s)):
    index = np.where(a.bl.s0 == a.spars_s[i])[0][0]
    xu_p = np.array([a.bu.g_p.x1_grid[index]])
    xl_p = np.array([a.bl.g_p.x1_grid[index]])
    plt.plot([xu_p, xl_p], [a.bu.g_p.x3(xu_p), a.bl.g_p.x3(xl_p)], 'b', lw=3)
    xu_c = np.array([a.bu.g.x1_grid[index]])
    xl_c = np.array([a.bl.g.x1_grid[index]])
    plt.plot([xu_c, xl_c], [a.bu.g.x3(xu_c), a.bl.g.x3(xl_c)], '.5', lw=3)

upper = np.loadtxt('case_study_7_upper.csv', delimiter=',')
lower = np.loadtxt('case_study_7_lower.csv', delimiter=',')
plt.scatter(upper[0, :], upper[1, :], c='.5', label='Abaqus', edgecolors='k',
            zorder=10, marker="^")
plt.scatter(lower[0, :], lower[1, :], c='.5', edgecolors='k', zorder=10,
            marker="^")
print('spars', a.bl.g.spar_x, a.bl.g.spar_y)
plt.scatter([a.bl.g.spar_x], [a.bl.g.spar_y], c='g', label='Lower spar', zorder=20, s=40)
# plt.axis('equal')

plt.figure()
dy = a.bu.g_p.x3(a.bu.g_p.x1_grid, diff='x1')
cos = dy/np.sqrt(1+dy*dy)
plt.plot(a.bu.g_p.x1_grid, dy, 'b',
         label='Upper Parent', lw=3)
plt.plot(a.bu.g_p.x1_grid, dy, '--b',
         label='Upper Parent', lw=3)
dy = a.bu.g.x3(a.bu.g.x1_grid, diff='x1')
cos = dy/np.sqrt(1+dy*dy)
plt.plot(a.bu.g.x1_grid, a.bu.g.x3(a.bu.g.x1_grid, diff='x1'), c='.5',
         label='Upper Child', lw=3)
plt.plot(a.bu.g.x1_grid, cos, '--', c='.5',
         label='Upper Child', lw=3)
dy = a.bl.g.x3(a.bl.g.x1_grid, diff='x1')
cos = dy/np.sqrt(1+dy*dy)
plt.plot(a.bl.g.x1_grid, a.bl.g.x3(a.bl.g.x1_grid, diff='x1'), c='g',
         label='Lower Child', lw=3)
plt.plot(a.bl.g.x1_grid, cos, '--', c='g',
         label='Lower Child', lw=3)

plt.scatter(upper[0, :], np.gradient(upper[1, :], upper[0, :]), c='.5', label='Abaqus Upper', edgecolors='k',
            zorder=10, marker="^")
plt.scatter(lower[0, :], np.gradient(lower[1, :], lower[0, :]), c='g', edgecolors='k', zorder=10,
            marker="s", label='Abaqus Lower')

plt.legend()
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
# plt.axis('equal')
# plt.legend()
# plt.gca().set_aspect('equal', adjustable='box')
plt.show()
