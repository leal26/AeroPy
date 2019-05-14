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

from aeropy.geometry.airfoil import CST
from aeropy.CST_2D import *


# Just as quick trick, to make upper morph I just mirror the image in regards to x
inverted = False
# Defines if basckwards or forwards morphing
morphing_direction = 'forwards'

# ==============================================================================
# Calculate dependent shape function parameters
# ==============================================================================


def calculate_dependent_shape_coefficients(Au_C_1_to_n,
                                           psi_spars, Au_P, Al_P, deltaz, c_P,
                                           morphing='backwards'):
    """Calculate  dependent shape coefficients for children configuration for a 4 order
    Bernstein polynomial and return the children upper, lower shape
    coefficients, children chord and spar thicknesses. _P denotes parent parameters"""
    def calculate_AC_u0(AC_u0):
        Au_C = [AC_u0] + Au_C_1_to_n
        c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz)
        return np.sqrt(c_P/c_C)*Au_P[0]

    # Bersntein Polynomial
    def K(r, n):
        K = math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K
    # Bernstein Polynomial order
    n = len(Au_C_1_to_n)

    # Find upper shape coefficient though iterative method since Au_0 is unknown
    # via fixed point iteration
    #AC_u0 = optimize.fixed_point(calculate_AC_u0, Au_P[0])
    # print AC_u0
    error = 9999
    AC_u0 = Au_P[0]
    while error > 1e-9:
        before = AC_u0
        AC_u0 = calculate_AC_u0(AC_u0)
        error = abs(AC_u0 - before)

    # Because the output is an array, need the extra [0]
    Au_C = [AC_u0] + Au_C_1_to_n

    # Now that AC_u0 is known we can calculate the actual chord and AC_l0
    c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz)
    AC_l0 = np.sqrt(c_P/c_C)*Al_P[0]
    # print(Au_C)
    # print(Au_P)
    # print(Al_P)
    # print(c_C, AC_l0, AC_u0)
    # print '0 lower shape coefficient: ',AC_l0
    # Calculate thicknessed and tensor B for the constraint linear system problem
    spar_thicknesses = []
    A0 = AC_u0 + AC_l0

    if morphing == 'backwards':
        b_list = np.zeros((n, 1))
        for j in range(len(psi_spars)):
            psi_j = psi_spars[j]
            # Calculate the spar thickness in meters from parent, afterwards, need to
            # adimensionalize for the goal airfoil by dividing by c_goal
            t_j = calculate_spar_distance(psi_spars[j], Au_C, Au_P, Al_P, deltaz, c_P)

            spar_thicknesses.append(t_j)
            b_list[j] = (t_j/c_C - psi_j*deltaz/c_C)/((psi_j**0.5)*(1-psi_j)) - A0*(1-psi_j)**n

        B = np.zeros((n, n))
        # j is the row dimension and i the column dimension in this case
        for j in range(n):
            for i in range(n):
                # Because in Python counting starts at 0, need to add 1 to be
                # coherent for equations
                r = i + 1
                B[j][i] = K(r, n)*(psi_spars[j]**r)*(1-psi_spars[j])**(n-r)

        A_bar = np.dot(inv(B), b_list)

        Al_C = [AC_l0]
        for i in range(len(A_bar)):
            Al_C.append(A_bar[i][0] - Au_C[i+1])  # extra [0] is necessary because of array

    elif morphing == 'forwards':
        f = np.zeros((n, 1))
        # psi/xi coordinates for lower surface of the children configuration
        psi_lower_children = []
        xi_lower_children = []
        xi_upper_children = []

        c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz)
        # print(c_C, AC_u0, AC_l0)
        # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
        psi_upper_children = []
        for j in range(len(psi_spars)):
            psi_upper_children.append(calculate_psi_goal(psi_spars[j], Au_P, Au_C, deltaz,
                                                         c_P, c_C))
        # Calculate xi for upper children. Do not care about lower so just gave it random shape coefficients
        xi_upper_children = CST(psi_upper_children, 1., deltasz=[
                                deltaz/2./c_C, deltaz/2./c_C],  Al=Au_C, Au=Au_C)
        xi_upper_children = xi_upper_children['u']

        # print xi_upper_children

        # Debugging section
        x = np.linspace(0, 1)
        y = CST(x, 1., deltasz=[deltaz/2./c_C, deltaz/2./c_C],  Al=Au_C, Au=Au_C)
        # plt.plot(x,y['u'])
        # plt.scatter(psi_upper_children, xi_upper_children)
        # plt.grid()
        # plt.show()
        # BREAK
        for j in range(len(psi_spars)):
            xi_parent = CST(psi_spars, 1., deltasz=[
                            deltaz/2./c_P, deltaz/2./c_P],  Al=Al_P, Au=Au_P)
            delta_j_P = xi_parent['u'][j]-xi_parent['l'][j]
            t_j = c_P*(delta_j_P)
            # Claculate orientation for children
            s_j = calculate_spar_direction(psi_spars[j], Au_P, Au_C, deltaz, c_C)
            psi_l_j = psi_upper_children[j]-delta_j_P/c_C*s_j[0]
            xi_l_j = xi_upper_children[j]-delta_j_P/c_C*s_j[1]

            spar_thicknesses.append(t_j)
            psi_lower_children.append(psi_l_j)
            xi_lower_children.append(xi_l_j)

            f[j] = (2*xi_l_j + psi_l_j*deltaz/c_C) / \
                (2*(psi_l_j**0.5)*(psi_l_j-1)) - AC_l0*(1-psi_l_j)**n

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

    Psi = np.array(x)/chord
    Xi = np.array(y)/chord
    EndThickness = EndThickness/chord

    # If A0 is given (constant leading edge assumption)
    if A0 is not None:
        n = len(x)
        T = np.zeros((n, n))
        t = np.zeros((n, 1))

        for j in range(1, n+1):
            jj = j - 1
            for i in range(1, n+1):
                ii = i - 1
                T[jj][ii] = K(i, n) * Psi[jj]**i * (1-Psi[jj])**(n-i)

            t[jj] = (Xi[jj] - Psi[jj]*EndThickness) / \
                (Psi[jj]**N1*(1-Psi[jj])**N2) - A0*(1-Psi[jj])**n
        # Calculate the inverse
        A = np.dot(inv(T), t)
        A = [A0] + list(A.transpose()[0])
    # If A0 is unknown
    else:
        n = len(x) - 1
        T = np.zeros((n+1, n+1))
        t = np.zeros((n+1, 1))
        for j in range(n+1):
            for i in range(n+1):
                T[j][i] = K(i, n) * Psi[j]**i * (1-Psi[j])**(n-i)

            t[j] = (Xi[j] - Psi[j]*EndThickness)/(Psi[j]**N1*(1-Psi[j])**N2)
        # Calculate the inverse
        A = np.dot(inv(T), t)
        A = list(A.transpose()[0])
    return A


def calculate_strains(Au_P, Al_P, c_P, Au_C, Al_C, c_C, deltaz, psi_spars, spar_thicknesses):
    # Calculate psi_flats (non-dimensional location of the itersection of
    # the spars with the lower surface
    psi_flats = []
    for j in range(len(psi_spars)):
        psi_parent_j = psi_spars[j]
        # Calculate psi at landing
        # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
        psi_children_j = calculate_psi_goal(psi_parent_j, Au_P, Au_C, deltaz, c_P, c_C)
        x_children_j = psi_children_j*c_C
        s = calculate_spar_direction(psi_spars[j], Au_P, Au_C, deltaz, c_C)
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


def plot_airfoil(AC, psi_spars, c_P, deltaz, Au_P, Al_P, image='plot',
                 iteration=0, return_coordinates=True, dir='current',
                 morphing_direction='backwards'):
    import matplotlib.pyplot as plt

    # plt.figure()
    n = len(Au_P) - 1
    Au_C, Al_C, c_C, spar_thicknesses = calculate_dependent_shape_coefficients(
        AC,
        psi_spars, Au_P, Al_P,
        deltaz, c_P, morphing=morphing_direction)
    print('CST chord', c_C)
    # ==============================================================================
    #  Plot results
    # ==============================================================================
    np.set_printoptions(precision=20)
    x = np.linspace(0, c_C, 1000)
    y = CST(x, c_C, deltasz=[deltaz/2., deltaz/2.],  Al=Al_C, Au=Au_C)
    plt.plot(x, y['u'], 'b', label='Children', lw=1)
    plt.plot(x, y['l'], '-b', label=None, lw=1)

    # store variables in case return_coordinates is True
    x = list(x[::-1]) + list(x[1:])
    y = list(y['u'][::-1]) + list(y['l'][1:])

    children_coordinates = {'x': x, 'y': y}
    x = np.linspace(0, c_P, 1000)
    y = CST(x, c_P, deltasz=[deltaz/2., deltaz/2.],  Al=Al_P, Au=Au_P)
    plt.plot(x, y['u'], 'r--', label='Parent', lw=1)
    plt.plot(x, y['l'], 'r--', label=None, lw=1)

    y_limits = y

    if morphing_direction == 'backwards':
        for i in range(len(psi_spars)):
            psi_i = psi_spars[i]
            # Calculate psi at landing
            psi_goal_i = calculate_psi_goal(psi_i, Au_C, Au_P, deltaz, c_C, c_P)
            x_goal_i = psi_goal_i*c_P
            # Calculate xi at landing
            temp = CST(x_goal_i, c_P, [deltaz/2., deltaz/2.], Al=Al_P, Au=Au_P)
            y_goal_i = temp['u']

            # calculate spar direction
            s = calculate_spar_direction(psi_i, Au_C, Au_P, deltaz, c_P)

            plt.plot([x_goal_i, x_goal_i - spar_thicknesses[i]*s[0]],
                     [y_goal_i, y_goal_i - spar_thicknesses[i]*s[1]], 'r--')

            y = CST(np.array([psi_i*c_C]), c_C, deltasz=[deltaz/2., deltaz/2.], Al=Al_C, Au=Au_C)
            plt.plot([psi_i*c_C, psi_i*c_C], [y['u'], y['u']-spar_thicknesses[i]], 'b', label=None)
    elif morphing_direction == 'forwards':
        for j in range(len(psi_spars)):
            psi_parent_j = psi_spars[j]
            # Calculate psi at landing
            # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
            psi_children_j = calculate_psi_goal(psi_parent_j, Au_P, Au_C, deltaz, c_P, c_C)
            x_children_j = psi_children_j*c_C

            # Calculate xi at landing
            temp = CST(x_children_j, c_C, [deltaz/2., deltaz/2.], Al=Al_C, Au=Au_C)
            y_children_j = temp['u']

            s = calculate_spar_direction(psi_spars[j], Au_P, Au_C, deltaz, c_C)

            # Print spars for children
            if not inverted:
                plt.plot([x_children_j, x_children_j - spar_thicknesses[j]*s[0]], [y_children_j,
                                                                                   y_children_j - spar_thicknesses[j]*s[1]], c='b', lw=1, label=None)
            else:
                plt.plot([x_children_j, x_children_j - spar_thicknesses[j]*s[0]],
                         [-y_children_j, -y_children_j + spar_thicknesses[j]*s[1]], c='b', lw=1, label=None)

            y = CST(np.array([psi_parent_j*c_P]), c_P,
                    deltasz=[deltaz/2., deltaz/2.], Al=Al_P, Au=Au_P)

            # Print spars for parents
            if not inverted:
                plt.plot([psi_parent_j*c_P, psi_parent_j*c_P],
                         [y['u'], y['u']-spar_thicknesses[j]], 'r--', lw=1, label=None)
            else:
                plt.plot([psi_parent_j*c_P, psi_parent_j*c_P], [-y['u'], -
                                                                y['u']+spar_thicknesses[j]], 'r--', lw=1, label=None)
                plt.plot([psi_i*c_C, psi_i*c_C], [y['u'], y['u'] -
                                                  spar_thicknesses[i]], 'b', label=None)
    plt.xlabel('$\psi$', fontsize=16)
    plt.ylabel(r'$\xi$', fontsize=16)
    plt.grid()
    plt.legend(loc="upper right")
    plt.gca().set_aspect(2, adjustable='box')
    x1, x2, y1, y2 = plt.axis()
    plt.axis((x1, x2, y1, 2*y2))

    # plt.axis([-0.005, c_L+0.005, min(y_limits['l'])-0.005, max(y_limits['l'])+0.01])
    if image == 'plot':
        plt.show()
    elif image == 'save':
        if dir == 'current':
            plt.savefig('%03i.pdf' % (iteration), bbox_inches='tight')
        else:
            cwd = os.getcwd()
            directory = os.path.join(cwd, dir)
            if not os.path.exists(directory):
                os.makedirs(directory)

            filename = os.path.join(directory, '%05i.png' % (iteration))
            plt.savefig(filename, bbox_inches='tight')
    if return_coordinates:
        return children_coordinates
