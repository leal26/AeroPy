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
from aeropy.CST_2D import calculate_c_baseline, calculate_psi_goal, calculate_spar_direction


def format_l(gu, gu_p, gl, gl_p):
    """Calculate  dependent shape coefficients for children configuration for a 4 order
    Bernstein polynomial and return the children upper, lower shape
    coefficients, children chord and spar thicknesses. _P denotes parent parameters"""
    def calculate_AC_u0(AC_u0):
        Au_C = [AC_u0] + Au_C_1_to_n
        c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz, N1=N1, N2=N2)
        return np.sqrt(c_P/c_C)*Au_P[0]

    psi_spars = [0.2, 0.3]
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
    deltaz = -gl.D[0]
    c_P = gu_p.chord
    # Find upper shape coefficient though iterative method since Au_0 is unknown
    # via fixed point iteration
    AC_u0 = gu.D[0]

    # Because the output is an array, need the extra [0]
    Au_C = [AC_u0] + Au_C_1_to_n

    # Now that AC_u0 is known we can calculate the actual chord and AC_l0
    c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz, N1=N1, N2=N2)
    AC_l0 = gl.D[0]
    # Calculate thicknessed and tensor B for the constraint linear system problem
    spar_thicknesses = []

    f = np.zeros((n, 1))
    # psi/xi coordinates for lower surface of the children configuration
    psi_lower_children = []
    xi_lower_children = []
    xi_upper_children = []

    c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz, N1=N1, N2=N2)
    psi_upper_children = []
    for j in range(len(psi_spars)):
        psi_upper_children.append(calculate_psi_goal(psi_spars[j], Au_P, Au_C, deltaz,
                                                     c_P, c_C, N1=N1, N2=N2))
    # Calculate xi for upper children. Do not care about lower so just gave it random shape coefficients
    xi_upper_children = CST(psi_upper_children, 1., deltasz=[
                            deltaz/c_C, deltaz/c_C], Al=Au_C, Au=Au_C, N1=N1, N2=N2)
    xi_upper_children = xi_upper_children['u']

    # print xi_upper_children

    for j in range(len(psi_spars)):
        xi_parent = CST(psi_spars, 1., deltasz=[
                        deltaz/c_P, deltaz/c_P], Al=Al_P, Au=Au_P, N1=N1, N2=N2)
        delta_j_P = xi_parent['u'][j]-xi_parent['l'][j]
        t_j = c_P*(delta_j_P)
        # Claculate orientation for children
        s_j = calculate_spar_direction(psi_spars[j], Au_P, Au_C, deltaz, c_C, N1=N1, N2=N2)
        psi_l_j = psi_upper_children[j]-delta_j_P/c_C*s_j[0]
        xi_l_j = xi_upper_children[j]-delta_j_P/c_C*s_j[1]

        spar_thicknesses.append(t_j)
        psi_lower_children.append(psi_l_j)
        xi_lower_children.append(xi_l_j)

        f[j] = (2*xi_l_j + psi_l_j*deltaz/c_C) / \
            (2*(psi_l_j**N1)*(psi_l_j-1)**N2) - AC_l0*(1-psi_l_j)**n

    F = np.zeros((n, n))
    # j is the row dimension and i the column dimension in this case
    for j in range(n):
        for i in range(n):
            # Because in Python counting starts at 0, need to add 1 to be
            # coherent for equations
            r = i + 1
            F[j][i] = K(r, n)*(psi_lower_children[j]**r)*(1-psi_lower_children[j])**(n-r)

    A_lower = np.dot(inv(F), f)

    Al_C = [AC_l0]
    for i in range(len(A_lower)):
        Al_C.append(A_lower[i][0])  # extra [0] is necessary because of array
    return Al_C + [-Al_C[0]]


def format_u(input, g=None, g_p=None):
    # COnsidering BC for zero derivative at the root
    return list(input) + [-input[0]]
    # return input


def format_input(input, gu=None, gu_p=None, gl=None, gl_p=None):
    Au = format_u(input[:3], gu, gu_p)
    gu_p.D[-1] = input[0]
    Al = format_u(input[3:], gl, gl_p)
    gl_p.D[-1] = input[3]
    print('Au', Au)
    print('Al', Al)
    # BREAK
    return Au, Al


g_upper = CoordinateSystem.CST(D=[0, 0, 0, 0], chord=1, color='b', N1=1, N2=1,
                               offset=.05)
g_lower = CoordinateSystem.CST(D=[0, 0, 0, 0], chord=1, color='b', N1=1, N2=1,
                               offset=-.05)
g_p = CoordinateSystem.CST(D=[0, 0, 0, 0], chord=1, color='k', N1=.5, N2=1)

s_upper = np.linspace(0, g_p.arclength(np.array([1]))[0], 51)
s_lower = np.linspace(0, g_p.arclength(np.array([1]))[0], 51)
p_upper = properties()
p_lower = properties()
l_upper = loads(concentrated_load=[[-np.sqrt(2)/2, -np.sqrt(2)/2]], load_s=[1])
l_lower = loads(concentrated_load=[[np.sqrt(2)/2, np.sqrt(2)/2]], load_s=[1-0.1])
a = coupled_beams(g_upper, g_lower, p_upper, p_lower, l_upper, l_lower, s_upper,
                  s_lower, ignore_ends=True)
a.calculate_x()

a.parameterized_solver(format_input=format_input, x0=list(g_upper.D[:-1]) + list(g_lower.D[:-1]))

# plt.figure()
# plt.plot(a.bu.g.x1_grid, a.bu.g.rho - a.bu.g_p.rho, label='Upper')
# plt.plot(a.bl.g.x1_grid, a.bl.g.rho - a.bl.g_p.rho, label='Lower')
# print('rho child', a.bl.g.rho)
# print('rho parent', a.bl.g_p.rho)
# plt.legend()

plt.figure()
plt.plot(a.bu.g_p.x1_grid, a.bu.g_p.x3(a.bu.g_p.x1_grid), 'b',
         label='Upper Parent', lw=3)
plt.plot(a.bu.g.x1_grid, a.bu.g.x3(a.bu.g.x1_grid), '.5',
         label='Upper Child: %.3f N' % -l_upper.concentrated_load[0][-1], lw=3)
plt.plot(a.bl.g_p.x1_grid, a.bl.g_p.x3(a.bl.g_p.x1_grid), 'b', linestyle='dashed',
         label='Lower Parent', lw=3)
plt.plot(a.bl.g.x1_grid, a.bl.g.x3(a.bl.g.x1_grid), '.5', linestyle='dashed',
         label='Lower Child: %.3f N' % -l_upper.concentrated_load[0][-1], lw=3)

# plt.scatter(abaqus_x, abaqus_y, c='.5', label='FEA: %.3f N' % -
#             l_upper.concentrated_load[0][-1], edgecolors='k', zorder=10, marker="^")
plt.legend()
plt.show()
