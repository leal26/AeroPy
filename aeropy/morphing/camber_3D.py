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
# from aeropy.CST_2D.module import CST_3D
from aeropy.CST_3D.module import CST_3D, calculate_c_baseline, \
    calculate_psi_goal, calculate_spar_direction

# TODO: Redo this for new 3D CST
# TODO: Make this object-oriented
# Just as quick trick, to make upper morph I just mirror in regards to x
inverted = False
# Defines if basckwards or forwards morphing
morphing_direction = 'forwards'

# ==============================================================================
# Calculate dependent shape function parameters
# ==============================================================================


def calculate_dependent_shape_coefficients(BP_p, BA_p, BP_c, chord_p, sweep_p,
                                           twist_p, delta_TE_p, sweep_c,
                                           twist_c, eta_sampling, psi_spars,
                                           morphing='camber'):
    """Calculate  dependent shape coefficients for children configuration for a
       4 order Bernstein polynomial and return the children upper, lower shape
       coefficients, children chord and spar thicknesses. _P denotes parent
       parameters"""
    def calculate_BP_c0(BP_c0, c_P, deltaz, j):
        Au_C = extract_A(BP_c, j)
        Au_P = extract_A(BP_p, j)
        c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz)
        BP_c[0][j] = np.sqrt(c_P/c_C)*Au_P[0]
        return BP_c[0][j], c_C

    def extract_A(B, j):
        # Extracting shape coefficient data for column j
        A = []
        for i in range(n+1):
            A.append(B[i][j])
        return A
    # Bersntein Polynomial

    def K(r, n):
        K = math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K

    # Bernstein Polynomial orders (n is for psi, and m for eta)
    n = len(BP_p) - 1
    m = len(BP_p[0]) - 1
    p = len(psi_spars)
    q = len(eta_sampling)
    # print p,q,n,m
    # Define chord, sweep, and twist functions for parent
    chord_p = CST(eta_sampling, chord_p['eta'][1], chord_p['initial'],
                  Au=chord_p['A'], N1=chord_p['N1'], N2=chord_p['N2'],
                  deltasLE=chord_p['final'])
    sweep_p = CST(eta_sampling, sweep_p['eta'][1], deltasz=sweep_p['final'],
                  Au=sweep_p['A'], N1=sweep_p['N1'], N2=sweep_p['N2'])
    chord_p = chord_p[::-1]
    sweep_p = sweep_p
    twist_p = CST(eta_sampling, twist_p['eta'][1], twist_p['initial'],
                  Au=twist_p['A'], N1=twist_p['N1'], N2=twist_p['N2'],
                  deltasLE=twist_p['final'])
    delta_TE_p = CST(eta_sampling, delta_TE_p['eta'][1], delta_TE_p['initial'],
                     Au=delta_TE_p['A'], N1=delta_TE_p['N1'],
                     N2=delta_TE_p['N2'], deltasLE=delta_TE_p['final'])
    # Initialize chord, sweep, and twist functions for child
    chord_c = []
    # Initialize child active matrix
    BA_c = []
    for i in range(n+1):
        temp = []
        for j in range(m+1):
            temp.append(0)
        BA_c.append(temp,)
    # Find upper shape coefficient though iterative method since Au_0 is
    # unknown via fixed point iteration
    for k in range(q):
        error = 9999
        BP_c0 = BP_p[0][k]
        while error > 1e-9:
            before = BP_c0
            c_P = chord_p[k]
            deltaz = delta_TE_p[k]
            [BP_c0, c_c] = calculate_BP_c0(BP_c0, c_P, deltaz, k)
            error = abs(BP_c0 - before)
        BP_c[0][k] = BP_c0
        BA_c[0][k] = np.sqrt(c_P/c_c)*BA_p[0][k]
        chord_c.append(c_c)
        print(c_c, BP_c0, np.sqrt(c_P/c_c)*BA_p[0][k])
    # Calculate thickness and tensor C for the constraint linear system problem
    psi_A_c = []

    if morphing == 'camber':
        f = np.zeros((q, p))
        for l in range(q):
            # Converting everything from 3D to 2D framework
            Au_P = extract_A(BP_p, l)
            Al_P = extract_A(BA_p, l)
            Au_C = extract_A(BP_c, l)
            c_P = chord_p[l]
            c_C = chord_c[l]
            deltaz = delta_TE_p[l]

            # psi/xi coordinates for lower surface of children configuration
            psi_lower_children = []
            xi_upper_children = []

            # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
            psi_upper_children = []
            for j in range(len(psi_spars)):
                psi_i = calculate_psi_goal(psi_spars[j], Au_P, Au_C, deltaz,
                                           c_P, c_C)
                psi_upper_children.append(psi_i)
            # Calculate xi for upper children. Do not care about lower so just
            # gave it random shape coefficients
            xi_upper_children = CST(psi_upper_children, 1., deltasz=[
                                    deltaz/2./c_C, deltaz/2./c_C], Al=Au_C,
                                    Au=Au_C)
            xi_upper_children = xi_upper_children['u']

            # print xi_upper_children

            # Debugging section
            # x = np.linspace(0, 1)
            # y = CST(x, 1., deltasz=[deltaz/2./c_C, deltaz/2./c_C],
            #         Al=Au_C, Au=Au_C)
            # plt.plot(x,y['u'])
            # plt.scatter(psi_upper_children, xi_upper_children)
            # plt.grid()
            # plt.show()
            # BREAK
            for k in range(len(psi_spars)):
                xi_parent = CST(psi_spars, 1., deltasz=[
                                deltaz/2./c_P, deltaz/2./c_P], Al=Al_P,
                                Au=Au_P)
                delta_k_P = xi_parent['u'][k]-xi_parent['l'][k]
                # t_k = c_P*(delta_k_P)
                # Claculate orientation for children
                s_k = calculate_spar_direction(psi_spars[k], Au_P, Au_C,
                                               deltaz, c_C)
                psi_l_k = psi_upper_children[k]-delta_k_P/c_C*s_k[0]
                xi_l_k = xi_upper_children[k]-delta_k_P/c_C*s_k[1]

                psi_lower_children.append(psi_l_k)

                f_y = 0
                for j in range(m+1):
                    f_y += BA_c[0][j]*(1-psi_l_k)**n \
                        * (K(j, m)*eta_sampling[l]**j
                           * (1-eta_sampling[l])**(m-j))

                f[l][k] = (2*xi_l_k + psi_l_k*deltaz/c_C) / \
                          (2*(psi_l_k**0.5) * (psi_l_k-1)) - f_y

                # print 'f_y',f_y, eta_sampling[l],m,j
            # Store new children psi values
            psi_A_c.append(psi_lower_children)

        # Initialize F (avoiding using numpy)
        F = np.zeros([q, p, m+1, n])
        # F = []
        # for l in range(q):
        # tempk = []
        # for k in range(p):
        # tempj = []
        # for j in range(m+1):
        # tempi = []
        # for i in range(n):
        # tempi.append(0.0)
        # tempj.append(tempi)
        # tempk.append(tempj)
        # F.append(tempk)

        # j is the row dimension and i the column dimension in this case
        for l in range(q):
            for k in range(p):
                for j in range(m+1):
                    for i in range(n):
                        # Because in Python counting starts at 0,
                        # need to add 1 to be coherent for equations
                        ii = i + 1
                        Sx = K(ii, n)*(psi_A_c[l][k]**ii) * \
                            (1-psi_A_c[l][k])**(n-ii)
                        Sy = K(j, m)*(eta_sampling[l]**j) * \
                            (1-eta_sampling[l])**(m-j)
                        F[l][k][j][i] = Sx*Sy

        # print len(F), len(F[0]), len(F[0][0]), len(F[0][0][0])

        # Unfolding tensor
        F_matrix = np.zeros((n**2, n**2))
        f_vector = np.zeros((n**2, 1))
        for l in range(q):
            for k in range(p):
                for j in range(m+1):
                    for i in range(n):
                        ii = n*(l) + k
                        jj = n*(j) + i

                        F_matrix[ii][jj] = F[l][k][j][i]
                        f_vector[ii] = f[l][k]

        solution = np.linalg.solve(F_matrix, f_vector)

        for j in range(m+1):
            for i in range(n):
                jj = n*(j) + i
                BA_c[i+1][j] = solution[jj][0]
        print(BA_c)
    return BA_c, chord_c
