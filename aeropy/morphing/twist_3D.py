# -*- coding: utf-8 -*-
"""
Objective: create an airfoil with a leading edge restriction, same upper length
restriction, othogonal upper spars and constant thicknesses in four places

Created on Mon Oct 17 10:36:34 2016

@author: Pedro
"""
from __future__ import print_function
import os
import math
import numpy as np
from numpy.linalg import inv

from aeropy.airfoil_module import CST
from aeropy.CST.module_2D import *

# Just as quick trick, to make upper morph I just mirror the image in regards to x
inverted = False
# Defines if basckwards or forwards morphing
morphing_direction = 'forwards'


def calculate_c_baseline(c_L, Au_C, Au_L, deltaz, l_LE=0, eps_LE=0, psi_P_u1=0):
    """Equations in the New_CST.pdf. Calculates the upper chord in order for
       the cruise and landing airfoils ot have the same length."""

    def integrand(psi, Au, delta_xi):
        return np.sqrt(1 + dxi_u(psi, Au, delta_xi)**2)

    def f(c_C):
        """Function dependent of c_C and that outputs c_C."""
        y_C, err = quad(integrand, 0, 1, args=(Au_C, deltaz/c_C))
        y_L, err = quad(integrand, psi_P_u1, 1, args=(Au_L, deltaz/c_L))
        y_LE, err = quad(integrand, 0, psi_P_u1, args=(Au_L, deltaz/c_L))
        return c_L*((1-eps_LE)*(l_LE+y_LE)+y_L)/y_C
    c_C = optimize.fixed_point(f, [c_L])
    # In case the calculated chord is really close to the original, but the
    # algorithm was not able to make them equal
    if abs(c_L - c_C) < 1e-7:
        return c_L
    # The output is an array so it needs the extra [0]
    return c_C[0]


def calculate_psi_goal(psi_baseline, Au_baseline, Au_goal, deltaz,
                       c_baseline, c_goal, l_LE, eps_LE, psi_1):
    """Find the value for psi that has the same location w on the upper
    surface of the goal as psi_baseline on the upper surface of the
    baseline"""

    def integrand(psi_baseline, Au, deltaz, c):
        return c*np.sqrt(1 + dxi_u(psi_baseline, Au, deltaz/c)**2)

    def equation(psi_goal, Au_goal, deltaz, c):
        if psi_goal != psi_1:
            L_baseline, err = quad(integrand, psi_1, psi_baseline, args=(Au_baseline, deltaz,
                                                                         c_baseline))
        else:
            L_baseline = 0
        L_LE, err = quad(integrand, 0, psi_1, args=(Au_baseline, deltaz,
                                                    c_baseline))
        y, err = quad(integrand, 0, psi_goal, args=(Au_goal, deltaz, c))
        return y - (1-eps_LE)*(L_LE+c_baseline*l_LE) - L_baseline

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y = fsolve(equation, psi_baseline, args=(Au_goal, deltaz,
                                                 c_goal))
    return y[0]


def calculate_A0_moving_LE(psi_baseline, psi_goal_0, Au_baseline, Au_goal, deltaz,
                           c_baseline, l_LE, eps_LE):
    """Find the value for A_P0^c that has the same arc length for the first bay
       as for the parent."""

    def integrand(psi_baseline, Al, deltaz, c):
        return c*np.sqrt(1 + dxi_u(psi_baseline, Al, deltaz/c)**2)

    def equation(A0, L_baseline, Au_goal, deltaz):
        Au_goal[0] = A0
        c = calculate_c_baseline(c_P, Au_goal, Au_baseline, deltaz/c_P, l_LE, eps_LE, psi_spars[0])
        y, err = quad(integrand, 0, psi_goal_0, args=(Au_goal, deltaz, c))
        print('y', y, y - (1-eps_LE)*L_baseline, A0, c)
        return y - (1-eps_LE)*(L_baseline - c*l_LE)
    L_baseline, err = quad(integrand, 0, psi_baseline[0], args=(Au_baseline, deltaz,
                                                                c_baseline))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y = fsolve(equation, Au_goal[0], args=(L_baseline, Au_goal, deltaz))
    return y[0]


def calculate_spar_direction(psi_baseline, Au_baseline, Au_goal, deltaz, c_goal, l_LE, eps_LE, psi_spars):
    """Calculate the direction of the spar component based on a location
    at the upper surface for the cruise airfoil."""
    # Calculate cruise chord
    c_baseline = calculate_c_baseline(c_goal, Au_baseline, Au_goal,
                                      deltaz, l_LE, eps_LE, psi_spars[0])
    # Calculate psi at goal arifoil
    psi_goal = calculate_psi_goal(psi_baseline, Au_baseline, Au_goal, deltaz,
                                  c_baseline, c_goal, l_LE, eps_LE, psi_spars[0])
    # non-normalized direction
    s = np.zeros(2)
    t = np.zeros(2)
#    t_norm = np.sqrt(1 + (dxi_u(psi_goal, Au_goal[0], Au_goal[1], deltaz))**2)

    cbeta = calculate_cbeta(psi_baseline, Au_baseline,
                            deltaz/c_baseline)
    sbeta = np.sqrt(1-cbeta**2)

    t[0] = 1
    t[1] = dxi_u(psi_goal, Au_goal, deltaz/c_goal)
    t_norm = np.sqrt(t[0]**2 + t[1]**2)
    t = (1./t_norm)*t
#    s[0] = t_norm*cbeta - dxi_u(psi_goal, Au_goal[0], Au_goal[1], deltaz)
#    s[1] =  1

    s[1] = t[1]*cbeta + t[0]*sbeta
    s[0] = (cbeta - s[1]*t[1])/t[0]
    return s

# ==============================================================================
# Calculate dependent shape function parameters
# ==============================================================================


def calculate_dependent_shape_coefficients(Au_C_1_to_n,
                                           psi_spars, Au_P, Al_P, deltaz, c_P,
                                           morphing='backwards', l_LE=0, eps_LE=0):
    """Calculate  dependent shape coefficients for children configuration for a 4 order
    Bernstein polynomial and return the children upper, lower shape
    coefficients, children chord and spar thicknesses. _P denotes parent parameters"""
    def calculate_AC_u0(AC_u0, constant_LE=True):
        Au_C = [AC_u0] + Au_C_1_to_n
        if constant_LE:
            return np.sqrt(c_P/c_C)*Au_P[0]
        else:
            return calculate_A0_moving_LE(psi_spars, psi_lower_children[0], Au_P, Au_C, deltaz,
                                          c_P, l_LE, eps_LE)
    # Bersntein Polynomial

    def K(r, n):
        K = math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K
    # Bernstein Polynomial order
    # In case of leading edge radius constraint
    n = len(Au_C_1_to_n)

    # Find upper shape coefficient though iterative method since Au_0 is unknown
    # via fixed point iteration
    #AC_u0 = optimize.fixed_point(calculate_AC_u0, Au_P[0])
    # print AC_u0
    error = 9999

    psi_lower_children = psi_spars
    Au_C = [Au_P[0]] + Au_C_1_to_n  # [Au_P[0]] +

    former_chord = c_P
    while error > 1e-5:
        former_Au_C = []
        for i in range(len(Au_C)):
            former_Au_C.append(Au_C[i])

        # Because the output is an array, need the extra [0]
        error_A0 = 999
        # Now that AC_u0 is known we can calculate the actual chord and AC_l0
        # c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz/c_P, l_LE, eps_LE, psi_spars[0])
        Au_C[0] = calculate_AC_u0(Au_C[0], constant_LE=False)
        Al_C0 = Au_C[0]
        c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz/c_P, l_LE, eps_LE, psi_spars[0])
        #Al_C0 = Au_C[0]
        # print '0 lower shape coefficient: ',AC_l0
        # Calculate thicknessed and tensor B for the constraint linear system problem
        spar_thicknesses = []

        if morphing == 'forwards':
            f = np.zeros((n, 1))
            # psi/xi coordinates for lower surface of the children configuration
            psi_lower_children = []
            xi_lower_children = []
            xi_upper_children = []

            # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
            psi_upper_children = []
            for j in range(len(psi_spars)):
                print(j)
                psi_upper_children.append(calculate_psi_goal(psi_spars[j], Au_P, Au_C, deltaz,
                                                             c_P, c_C, l_LE, eps_LE, psi_spars[0]))

            # Calculate xi for upper children. Do not care about lower so just gave it random shape coefficients
            xi_upper_children = CST(psi_upper_children, 1., deltasz=[
                                    deltaz/2./c_C, deltaz/2./c_C],  Al=Au_C, Au=Au_C)
            xi_upper_children = xi_upper_children['u']

            # print xi_upper_children

            # Debugging section
            # x = np.linspace(0,1)
            # y = CST(x, 1., deltasz= [deltaz/2./c_C, deltaz/2./c_C],  Al= Au_C, Au =Au_C)
            # plt.plot(x,y['u'])
            # plt.scatter(psi_upper_children, xi_upper_children)
            # plt.grid()
            # plt.show()
            # BREAK
            print(Au_P, Au_C, len(psi_spars), n)
            for j in range(len(psi_spars)):
                xi_parent = CST(psi_spars, 1., deltasz=[
                                deltaz/2./c_P, deltaz/2./c_P],  Al=Al_P, Au=Au_P)
                delta_j_P = xi_parent['u'][j]-xi_parent['l'][j]
                t_j = c_P*(delta_j_P)
                # Claculate orientation for children
                s_j = calculate_spar_direction(
                    psi_spars[j], Au_P, Au_C, deltaz, c_C, l_LE, eps_LE, psi_spars)
                psi_l_j = psi_upper_children[j]-delta_j_P/c_C*s_j[0]
                xi_l_j = xi_upper_children[j]-delta_j_P/c_C*s_j[1]

                spar_thicknesses.append(t_j)
                psi_lower_children.append(psi_l_j)
                xi_lower_children.append(xi_l_j)

                f[j] = (2*xi_l_j + psi_l_j*deltaz/c_C) / \
                    (2*(psi_l_j**0.5)*(psi_l_j-1)) - Al_C0*(1-psi_l_j)**n

            F = np.zeros((n, n))
            # j is the row dimension and i the column dimension in this case
            for j in range(n):
                for i in range(n):
                    # Because in Python counting starts at 0, need to add 1 to be
                    # coherent for equations
                    r = i + 1
                    F[j][i] = K(r, n)*(psi_lower_children[j]**r)*(1-psi_lower_children[j])**(n-r)
            print(F)
            print(f)
            A_lower = np.dot(inv(F), f)
            print('result', A_lower)
            Al_C = [Al_C0]
            for i in range(len(A_lower)):
                Al_C.append(A_lower[i][0])  # extra [0] is necessary because of array
        error_denominator = 0
        print('before', former_Au_C, Au_C)
        for i in range(len(Au_C)):
            error_denominator += Au_C[i]**2
        error = 0
        for i in range(len(Al_C)):
            error += (former_Au_C[i] - Au_C[i])**2/error_denominator
        error = math.sqrt(error)
        # error = abs(c_C-former_chord)/c_C
        # AC_u0 = calculate_AC_u0(AC_u0, constant_LE=False)
        print(error, Al_C, Au_C)
        # former_chord = c_C
    return Au_C, Al_C, c_C, spar_thicknesses


def calculate_shape_coefficients_tracing(A0, x, y, N1, N2, chord=1., EndThickness=0):
    """
    inputs:
        - tip_displacement: {'x': value, 'y': value}
        - other_points: {'x': value, 'y': value}
        - A0: float value for first shape coefficient. Usually related to a constraint.
    """
    # Bersntein Polynomial
    def K(r, n):
        K = math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K

    n = len(x)

    print(x)
    Psi = np.array(x)/chord
    Xi = np.array(y)/chord

    EndThickness = EndThickness/chord
    T = np.zeros((n, n))
    t = np.zeros((n, 1))
    for j in range(1, n+1):
        jj = j - 1
        for i in range(1, n+1):
            ii = i - 1
            T[jj][ii] = K(i, n) * Psi[jj]**i * (1-Psi[jj])**(n-i)
        print(Xi[jj], EndThickness, Psi[jj], A0, Psi[jj]**N1*(1-Psi[jj])**N2)
        t[jj] = (Xi[jj] - Psi[jj]*EndThickness)/(Psi[jj]**N1*(1-Psi[jj])**N2) - A0*(1-Psi[jj])**n
    # Calculate the inverse
    A = np.dot(inv(T), t)
    A = [A0] + list(A.transpose()[0])
    return A


def calculate_strains(Au_P, Al_P, c_P, Au_C, Al_C, c_C, deltaz, psi_spars, spar_thicknesses):
    # Calculate psi_flats (non-dimensional location of the itersection of
    # the spars with the lower surface
    psi_flats = []
    for j in range(len(psi_spars)):
        psi_parent_j = psi_spars[j]
        # Calculate psi at landing
        # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
        psi_children_j = calculate_psi_goal(
            psi_parent_j, Au_P, Au_C, deltaz, c_P, c_C, l_LE, eps_LE, psi_spars[0])
        x_children_j = psi_children_j*c_C
        s = calculate_spar_direction(psi_spars[j], Au_P, Au_C, deltaz, c_C, l_LE, eps_LE, psi_spars)
        psi_flats.append(x_children_j - spar_thicknesses[j]*s[0])

    # Calculate initial lengths
    initial_lengths = []
    psi_list = [0.] + psi_spars + [c_P]
    for i in range(len(psi_list)-1):
        initial_lengths.append(calculate_arc_length(psi_list[i], psi_list[i+1], Al_P, deltaz, c_P))

    # Calculate final lengths
    final_lengths = []
    psi_list = [0.] + psi_flats + [c_C]  # In P configuration
    for i in range(len(psi_list)-1):
        final_lengths.append(calculate_arc_length(
            psi_list[i]*c_P/c_C, psi_list[i+1]*c_P/c_C, Al_C, deltaz, c_C))

    # Calculate strains
    strains = []

    for i in range(len(final_lengths)):
        strains.append((final_lengths[i]-initial_lengths[i])/initial_lengths[i])
    av_strain = (sum(final_lengths)-sum(initial_lengths))/sum(initial_lengths)
    # for i in range(len(strains)):
    # print 'Initial length: ' + str(initial_lengths[i]) + ', final length: ' + str(final_lengths[i]) + ', strains: ' + str(strains[i])

    return strains, av_strain


def plot_airfoil(AC, psi_spars, c_L, deltaz, Au_L, Al_L, image='plot',
                 iteration=0, return_coordinates=True, dir='current'):
    import matplotlib.pyplot as plt

    plt.figure()
    n = len(Au_L) - 1
    Au_C, Al_C, c_C, spar_thicknesses = calculate_dependent_shape_coefficients(
        AC,
        psi_spars, Au_L, Al_L,
        deltaz, c_L, morphing=morphing_direction)

    # ==============================================================================
    #  Plot results
    # ==============================================================================
    np.set_printoptions(precision=20)
    x = np.linspace(0, c_C, 1000)
    y = CST(x, c_C, deltasz=[deltaz/2., deltaz/2.],  Al=Al_C, Au=Au_C)
    plt.plot(x, y['u'], 'b', label='Children')
    plt.plot(x, y['l'], '-b', label=None)

    # store variables in case return_coordinates is True
    x = list(x[::-1]) + list(x[1:])
    y = list(y['u'][::-1]) + list(y['l'][1:])

    children_coordinates = {'x': x, 'y': y}
    x = np.linspace(0, c_L, 1000)
    y = CST(x, c_L, deltasz=[deltaz/2., deltaz/2.],  Al=Al_L, Au=Au_L)
    plt.plot(x, y['u'], 'r--', label='Parent')
    plt.plot(x, y['l'], 'r--', label=None)

    y_limits = y

    for i in range(len(psi_spars)):
        psi_i = psi_spars[i]
        # Calculate psi at landing
        psi_goal_i = calculate_psi_goal(psi_i, Au_C, Au_L, deltaz, c_C, c_L)
        x_goal_i = psi_goal_i*c_L
        # Calculate xi at landing
        temp = CST(x_goal_i, c_L, [deltaz/2., deltaz/2.], Al=Al_L, Au=Au_L)
        y_goal_i = temp['u']

        # calculate spar direction
        s = calculate_spar_direction(psi_i, Au_C, Au_L, deltaz, c_L)

        plt.plot([x_goal_i, x_goal_i - spar_thicknesses[i]*s[0]],
                 [y_goal_i, y_goal_i - spar_thicknesses[i]*s[1]], 'r--')

        y = CST(np.array([psi_i*c_C]), c_C, deltasz=[deltaz/2., deltaz/2.], Al=Al_C, Au=Au_C)
        plt.plot([psi_i*c_C, psi_i*c_C], [y['u'], y['u']-spar_thicknesses[i]], 'b', label=None)

    plt.xlabel('$\psi$', fontsize=16)
    plt.ylabel(r'$\xi$', fontsize=16)
    plt.grid()
    plt.legend(loc="upper right")
    plt.gca().set_aspect('equal', adjustable='box')
    x1, x2, y1, y2 = plt.axis()
    plt.axis((x1, x2, y1, 2*y2))

    # plt.axis([-0.005, c_L+0.005, min(y_limits['l'])-0.005, max(y_limits['l'])+0.01])
    if image == 'plot':
        plt.show()
    elif image == 'save':
        if dir == 'current':
            plt.savefig('%03i.png' % (iteration), bbox_inches='tight')
        else:
            cwd = os.getcwd()
            directory = os.path.join(cwd, dir)
            if not os.path.exists(directory):
                os.makedirs(directory)

            filename = os.path.join(directory, '%05i.png' % (iteration))
            plt.savefig(filename, bbox_inches='tight')
    if return_coordinates:
        return children_coordinates