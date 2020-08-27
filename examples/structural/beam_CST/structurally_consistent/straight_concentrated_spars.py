import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
import math
import numpy as np
from numpy.linalg import inv

from aeropy.structural.beam import beam_chen, coupled_beams
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem
from aeropy.geometry.airfoil import CST
from aeropy.CST_2D import calculate_c_baseline, calculate_psi_goal, calculate_spar_direction, S


def format_l(gu, gu_p, gl, gl_p):
    """Calculate  dependent shape coefficients for children configuration for a 4 order
    Bernstein polynomial and return the children upper, lower shape
    coefficients, children chord and spar thicknesses. _P denotes parent parameters"""
    def calculate_AC_u0(AC_u0):
        Au_C = [AC_u0] + Au_C_1_to_n
        c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz, N1=N1, N2=N2)
        return np.sqrt(c_P/c_C)*Au_P[0]

    # Bersntein Polynomial

    def K(r, n):
        K = math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K
    # Bernstein Polynomial order
    n = len(gl.D) - 2
    N1 = gl.N1
    N2 = gl.N2
    Au_C_1_to_n = gu.D[1:-1]
    Au_P = gu_p.D[:-1]
    Al_P = gl_p.D[:-1]
    deltaz_C = -gl.D[0]*gl.chord
    deltaz_P = -gl_p.D[0]*gl_p.chord
    c_P = gu_p.chord
    # Find upper shape coefficient though iterative method since Au_0 is unknown
    # via fixed point iteration
    AC_u0 = gu.D[0]

    # Because the output is an array, need the extra [0]
    Au_C = [AC_u0] + Au_C_1_to_n

    # Now that AC_u0 is known we can calculate the actual chord and AC_l0
    c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz_P, deltaz_C, N1=N1, N2=N2)

    Al_C = np.zeros(n+1)
    Al_C[0] = gl.D[0]
    Al_C[1+len(psi_spars):] = gl.D[1+len(psi_spars):-1]
    # print('test', gl.D[1+len(psi_spars):-1])
    # Calculate thicknessed and tensor B for the constraint linear system problem
    spar_thicknesses = []

    f = np.zeros((m, 1))
    # psi/xi coordinates for lower surface of the children configuration
    psi_lower_children = []
    xi_lower_children = []
    xi_upper_children = []

    c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz_P, deltaz_C, N1=N1, N2=N2)
    psi_upper_children = []
    for j in range(len(psi_spars)):
        psi_upper_children.append(calculate_psi_goal(psi_spars[j], Au_P, Au_C, deltaz_C,
                                                     c_P, c_C, N1=N1, N2=N2, ))
    # Calculate xi for upper children. Do not care about lower so just gave it random shape coefficients
    xi_upper_children = CST(psi_upper_children, 1., deltasz=[
                            deltaz_C/c_C, deltaz_C/c_C], Al=Au_C, Au=Au_C, N1=N1, N2=N2)
    xi_upper_children = xi_upper_children['u']

    # print xi_upper_children
    gl.spar_directions = []

    for j in range(len(psi_spars)):
        xi_parent = CST(psi_spars, 1., deltasz=[
                        deltaz_P/c_P, deltaz_P/c_P], Al=Al_P, Au=Au_P, N1=N1, N2=N2)
        delta_j_P = xi_parent['u'][j]-xi_parent['l'][j] + gu_p.offset/c_P - gl_p.offset/c_P

        t_j = c_P*(delta_j_P)
        # Claculate orientation for children
        s_j = calculate_spar_direction(
            psi_spars[j], Au_P, Au_C, deltaz_P, c_C, N1=N1, N2=N2, deltaz_goal=deltaz_C)
        psi_l_j = psi_upper_children[j]-delta_j_P/c_C*s_j[0]
        xi_l_j = xi_upper_children[j]+gu.offset-delta_j_P/c_C*s_j[1]

        spar_thicknesses.append(t_j)
        psi_lower_children.append(psi_l_j)
        xi_lower_children.append(xi_l_j)

        f[j] = (xi_l_j - psi_l_j*deltaz_C/c_C - gl.offset/c_C) / \
            ((psi_l_j**N1)*(1-psi_l_j)**N2) - Al_C[0]*(1-psi_l_j)**n

        for k in range(1+len(psi_spars), n):
            f[j] -= Al_C[k]*S(k, n, psi_l_j)
        gl.spar_directions.append(s_j)
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


upper = np.loadtxt('upper_beam_spar2.csv', delimiter=',')
lower = np.loadtxt('lower_beam_spar2.csv', delimiter=',')

psi_spars = [0.2]
m = len(psi_spars)

g_upper = CoordinateSystem.CST(D=[0., 0., 0., 0., 0., 0., 0., ], chord=1, color='b', N1=1, N2=1,
                               offset=.05)
g_lower = CoordinateSystem.CST(D=[0., 0., 0., 0., 0., 0., 0., ], chord=1, color='b', N1=1, N2=1,
                               offset=-.05)
g_p = CoordinateSystem.CST(D=[0., 0., 0., 0., 0., 0., 0., ], chord=1, color='k', N1=1, N2=1)

s_upper = np.linspace(0, g_p.arclength(np.array([1]))[0], 51)
s_lower = np.linspace(0, g_p.arclength(np.array([1]))[0], 51)
p_upper = properties()
p_lower = properties()
l_upper = loads(concentrated_load=[[-np.sqrt(2)/2, -np.sqrt(2)/2]], load_s=[1])
l_lower = loads(concentrated_load=[[np.sqrt(2)/2, np.sqrt(2)/2]], load_s=[1-0.1])
a = coupled_beams(g_upper, g_lower, p_upper, p_lower, l_upper, l_lower, s_upper,
                  s_lower, ignore_ends=True, spars_s=psi_spars)
a.calculate_x()

a.parameterized_solver(format_input=format_input, x0=list(
    g_upper.D[:-1]) + list(g_lower.D[:1]) + list(g_lower.D[2:-1]))
print('LOADS', a.bl.l.concentrated_load, a.bu.l.concentrated_load)
plt.figure()
plt.plot(a.bu.g.x1_grid, a.bu.M, label='Upper')
plt.plot(a.bl.g.x1_grid, a.bl.M, label='Lower')
plt.legend()

index = np.where(a.bl.s == a.spars_s[0])[0][0]
plt.figure()
plt.plot(a.bu.g_p.x1_grid, a.bu.g_p.x3(a.bu.g_p.x1_grid), 'b',
         label='Upper Parent', lw=3)
plt.plot(a.bu.g.x1_grid, a.bu.g.x3(a.bu.g.x1_grid), '.5',
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
plt.scatter(a.bl.g.spar_psi, a.bl.g.spar_xi, c='g', label='spars')
plt.legend()
plt.show()
