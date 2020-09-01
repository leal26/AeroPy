import matplotlib.pyplot as plt
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


def format_l(gu, gu_p, gl, gl_p):
    """Calculate  dependent shape coefficients for children configuration for a 4 order
    Bernstein polynomial and return the children upper, lower shape
    coefficients, children chord and spar thicknesses. _P denotes parent parameters"""

    # Bersntein Polynomial

    def K(r, n):
        K = math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K
    # Bernstein Polynomial order
    n = len(gl.D) - 2
    N1 = gl.N1
    N2 = gl.N2
    Au_C = gu.D[0:-1]
    Au_P = gu_p.D[:-1]
    Al_P = gl_p.D[:-1]
    zl_P = -gl_p.D[0]*gl_p.chord
    zu_P = -gu_p.D[0]*gl_p.chord
    cu_P = gu_p.chord
    cl_P = gl_p.chord

    Al_C = np.zeros(n+1)
    Al_C[0] = gl.D[0]
    Al_C[1+len(psi_spars):] = gl.D[1+len(psi_spars):-1]
    # print('test', gl.D[1+len(psi_spars):-1])
    # Calculate thicknessed and tensor B for the constraint linear system problem
    f = np.zeros((m, 1))
    # psi/xi coordinates for lower surface of the children configuration

    xi_upper_children = []
    gu.internal_variables(arc_upper)
    cu_C = gu.chord
    # cu_C = calculate_c_baseline(cu_P, Au_P, Au_C, zu_P, zu_C, N1=N1, N2=N2)
    zu_C = -gu.D[0]*cu_C

    psi_upper_children = []
    for j in range(len(psi_spars)):
        psi_upper_children.append(calculate_psi_goal(psi_spars[j], Au_P, Au_C, zu_P,
                                                     cu_P, cu_C, N1=N1, N2=N2, deltaz_goal=zu_C))
    # Calculate xi for upper children. Do not care about lower so just gave it random shape coefficients
    xi_upper_children = CST(psi_upper_children, 1., deltasz=[
        zu_C/cu_C, -zu_C/cu_C], Al=Au_C, Au=Au_C, N1=N1, N2=N2)
    xi_upper_children = xi_upper_children['u']

    gl.spar_psi_upper = psi_upper_children
    gl.spar_xi_upper = np.array(xi_upper_children) + gu.offset/cu_C
    # print xi_upper_children
    gl.spar_directions = []
    gl.delta_P = []
    for j in range(len(psi_spars)):
        xi_parent = CST(psi_spars, 1., deltasz=[
            zu_P/cu_P, -zl_P/cl_P], Al=Al_P, Au=Au_P, N1=N1, N2=N2)
        delta_j_P = cu_P*(xi_parent['u'][j]-xi_parent['l'][j]) + gu_p.offset - gl_p.offset

        # Claculate orientation for children
        s_j = calculate_spar_direction(
            psi_spars[j], Au_P, Au_C, zu_P, cu_C, N1=N1, N2=N2, deltaz_goal=zu_C)
        gl.spar_directions.append(s_j)
        gl.delta_P.append(delta_j_P)
    # print('spar_directions', gl.spar_directions)
    error = 999
    gl.internal_variables(arc_lower)
    chord0 = gl.chord
    cl_C = gl.chord
    while error > 1e-6:
        # cl_C = calculate_c_baseline(cl_P, Al_P, Al_C, zl_P, zl_C, N1=N1, N2=N2)
        zl_C = -gl.D[0]*cl_C
        psi_lower_children = []
        xi_lower_children = []

        for j in range(len(psi_spars)):
            delta_j_P = gl.delta_P[j]
            s_j = gl.spar_directions[j]
            psi_l_j = (psi_upper_children[j]*cu_C-delta_j_P*s_j[0])/cl_C
            xi_l_j = (xi_upper_children[j]*cu_C+gu.offset-delta_j_P*s_j[1])/cl_C
            # print('INSIDE', xi_upper_children[j], cu_C, gu.offset, -delta_j_P*s_j[1], cl_C)
            psi_lower_children.append(psi_l_j)
            xi_lower_children.append(xi_l_j)

            f[j] = (xi_l_j - psi_l_j*zl_C/cl_C - gl.offset/cl_C) / \
                ((psi_l_j**N1)*(1-psi_l_j)**N2) - Al_C[0]*(1-psi_l_j)**n

            for k in range(1+len(psi_spars), n):
                f[j] -= Al_C[k]*S(k, n, psi_l_j)

        gl.spar_psi = psi_lower_children
        gl.spar_xi = xi_lower_children
        F = np.zeros((m, m))
        # j is the row dimension and i the column dimension in this case
        for j in range(m):
            for i in range(m):
                # Because in Python counting starts at 0, need to add 1 to be
                # coherent for equations
                r = i + 1
                F[j][i] = K(r, n)*(psi_lower_children[j]**r)*(1-psi_lower_children[j])**(n-r)

        A_lower = np.dot(inv(F), f)

        for i in range(len(A_lower)):
            Al_C[i+1] = A_lower[i][0]  # extra [0] is necessary because of array
        gl.D = list(Al_C) + [-Al_C[0]]
        gl.internal_variables(arc_lower)
        cl_C = gl.chord
        error = abs(chord0-cl_C)
        chord0 = cl_C
        # print('A in', list(Al_C) + [-Al_C[0]])
        # print('c', cl_C, cu_C, error)
        # print('error', error)
    # print('directions', gl.spar_directions)
    # print('coord', gl.spar_psi, gl.spar_xi)
    # print('Coefficients', Al_C)
    # BREAK
    return list(Al_C) + [-Al_C[0]]


def format_u(input, g=None, g_p=None):
    # COnsidering BC for zero derivative at the root
    return list(input) + [-input[0]]
    # return input


def format_input(input, gu=None, gu_p=None, gl=None, gl_p=None):
    Au = format_u(input[:len(gu.D)-1], gu, gu_p)
    mn = len(gu.D) - 2 - m

    gu.D = Au
    gl.D[0] = input[len(gu.D)-1]
    gl.D[-mn-1:-1] = input[-mn:]
    gl.D[-1] = -input[len(gu.D)-1]
    # print('D', gl.D, input[-mn:], mn, input)
    # BREAK
    #
    # print('mn')
    Al = format_l(gu, gu_p, gl, gl_p)
    # gl_p.D[-1] = input[3]
    # print('input', input)
    # print('Au', Au)
    # print('Al', Al)

    # BREAK
    return Au, Al


def constraint_f(input):
    if hasattr(a.bl.g, 'spar_psi'):
        return arc_spar - a.bl.g.arclength(np.array(a.bl.g.spar_psi))[0]
    else:
        return 0


warnings.filterwarnings("ignore", category=RuntimeWarning)

upper = np.loadtxt('upper_beam_spar2.csv', delimiter=',')
lower = np.loadtxt('lower_beam_spar2.csv', delimiter=',')

psi_spars = [0.25]
m = len(psi_spars)

g_upper = CoordinateSystem.CST(D=[0., 0., 0., 0., 0., 0., 0., 0.], chord=1, color='b', N1=1, N2=1,
                               offset=.05)
g_lower = CoordinateSystem.CST(D=[0., 0., 0., 0., 0., 0., 0., 0.], chord=1, color='b', N1=1, N2=1,
                               offset=-.05)
g_p = CoordinateSystem.CST(D=[0., 0., 0., 0., 0., 0., 0., 0.], chord=1, color='k', N1=1, N2=1)

s_upper = np.linspace(0, g_p.arclength(np.array([1]))[0], 21)
s_lower = np.linspace(0, g_p.arclength(np.array([1]))[0], 21)
p_upper = properties()
p_lower = properties()
l_upper = loads(concentrated_load=[[-np.sqrt(2)/2, -np.sqrt(2)/2]], load_s=[1])
l_lower = loads(concentrated_load=[[np.sqrt(2)/2, np.sqrt(2)/2]], load_s=[1-0.1])
# l_upper = loads(concentrated_load=[[0, -1]], load_s=[1])
# l_lower = loads(concentrated_load=[[0, 1]], load_s=[1])
arc_upper = 1.0
arc_lower = 1.0
arc_spar = g_lower.arclength(np.array([psi_spars[0]]))[0]

a = coupled_beams(g_upper, g_lower, p_upper, p_lower, l_upper, l_lower, s_upper,
                  s_lower, ignore_ends=True, spars_s=psi_spars)
a.calculate_x()
constraints = ({'type': 'eq', 'fun': constraint_f})
# a.formatted_residual(format_input=format_input, x0=[
#                      .1, .1, .1, -.1] + list(g_lower.D[2:-1]))
# a.formatted_residual(format_input=format_input, x0=list(
# g_upper.D[:-1]) + list(g_lower.D[:1]) + list(g_lower.D[2:-1]))
a.parameterized_solver(format_input=format_input, x0=list(
    g_upper.D[:-1]) + list(g_lower.D[:-1]))
print('loads', a.bl.l.concentrated_load, a.bu.l.concentrated_load)
plt.figure()
plt.plot(a.bu.g.x1_grid, a.bu.M, label='Upper')
plt.plot(a.bl.g.x1_grid, a.bl.M, label='Lower')
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
x = [a.bu.g.chord*a.bl.g.spar_psi_upper[0], a.bl.g.chord*a.bl.g.spar_psi[0]]
y = [a.bu.g.chord*a.bl.g.spar_xi_upper[0], a.bl.g.chord*a.bl.g.spar_xi[0]]
dx = x[1]-x[0]
dy = y[1]-y[0]
norm = math.sqrt(dx**2+dy**2)
# print('spar direction', a.bl.g.spar_directions)
# print('actual direction', dx/norm, dy/norm)
plt.plot(x, y, c='g', label='spars', lw=3)
plt.arrow(x[0], y[0], -a.bl.g.spar_directions[0][0]*a.bl.g.delta_P[0],
          -a.bl.g.spar_directions[0][1]*a.bl.g.delta_P[0])
# print(a.bl.g.delta_P[0])
plt.legend()
# plt.gca().set_aspect('equal', adjustable='box')
plt.show()
