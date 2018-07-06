'''Programmed by Antoine BALDO'''
from __future__ import print_function

import sys
import numpy as np
import matplotlib.pyplot as plt

from aeropy.morphing.camber_2D import calculate_shape_coefficients_tracing, calculate_dependent_shape_coefficients
from aeropy.geometry.airfoil import CST
from aeropy.CST.module_2D import *

morphing_direction = 'forwards'
inverted = False

# Points where children should be
ValX = [0.1 ,0.25,0.42,0.65,.9]
ValY = [0.02,0.03,0.04,0.05,.05]

# Geoemtric propeties for parent
c_P = 1.
Au_P =  [0.10887, 0.1187, 0.07843, 0.12084, 0.07919, 0.09840]
Al_P =  [0.11117, 0.1000, 0.1239, 0.06334, 0.11539, 0.10400]
psi_spars = [0.2,0.3,0.5,0.7,0.9]
deltaz = 0

# Initialize values before iteration
AC_u0 = Au_P[0]
c_C = c_P
# Iterative method is necessary because the value of Au_C is not known
tol = 1e-6
error = 99999.
counter = 1
while error>tol:
    # tracing :
    A = calculate_shape_coefficients_tracing(AC_u0, ValX, ValY, 0.5, 1.,c_C, deltaz)
    # structurally_consistent :
    Au_C, Al_C, c_C, spar_thicknesses = calculate_dependent_shape_coefficients(
                                                        A[1:], psi_spars, Au_P, Al_P,
                                                        deltaz, c_P, morphing=morphing_direction)
    error = abs((AC_u0-Au_C[0])/AC_u0) 
    print('Iteration: ' + str(counter) + ', Error: ' +str(error))
    AC_u0 = Au_C[0]
    counter += 1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x = np.linspace(0, c_C, 1000)
y = CST(x, c_C, deltasz= [deltaz/2., deltaz/2.],  Al= Al_C, Au =Au_C)
plt.plot(x, y['u'],'b',label = 'Children', lw=2)
plt.plot(x, y['l'],'b',label = None, lw=2)

# Print shape for parent
x = np.linspace(0, c_P, 1000)
y = CST(x, c_P, deltasz= [deltaz/2., deltaz/2.],  Al= Al_P, Au =Au_P)
plt.plot(x, y['u'],'r--',label='Parent', lw=2)
plt.plot(x, y['l'],'r--',label = None, lw=2)

if morphing_direction == 'forwards':
    psi_flats = []
    intersections_x_children = [0]
    intersections_y_children = [0]
    intersections_x_parent = [0]
    intersections_y_parent = [0]
    for j in range(len(psi_spars)):
        psi_parent_j = psi_spars[j]
        # Calculate psi at landing
        # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
        psi_children_j = calculate_psi_goal(psi_parent_j, Au_P, Au_C, deltaz, c_P, c_C)
        x_children_j = psi_children_j*c_C

        # Calculate xi at landing
        temp = CST(x_children_j, c_C, [deltaz/2., deltaz/2.], Al= Al_C, Au =Au_C)
        y_children_j = temp['u']

        s = calculate_spar_direction(psi_spars[j], Au_P, Au_C, deltaz, c_C)

        # Print spars for children
        plt.plot([x_children_j, x_children_j - spar_thicknesses[j]*s[0]],[y_children_j, y_children_j - spar_thicknesses[j]*s[1]], c = 'b', lw=2, label=None)
        psi_flats.append(x_children_j - spar_thicknesses[j]*s[0])
        y = CST(np.array([psi_parent_j*c_P]), c_P, deltasz=[deltaz/2., deltaz/2.], Al= Al_P, Au =Au_P)

        intersections_x_children.append(x_children_j - spar_thicknesses[j]*s[0])
        intersections_y_children.append(y_children_j - spar_thicknesses[j]*s[1])

        # Print spars for parents
        plt.plot([psi_parent_j*c_P, psi_parent_j*c_P], [y['u'], y['u']-spar_thicknesses[j]], 'r--', lw=2, label = None)

        intersections_x_parent.append(psi_parent_j*c_P)
        intersections_y_parent.append(y['u']-spar_thicknesses[j])

plt.scatter([0]+ValX, [0]+ValY)
plt.xlabel('$\psi^p$', fontsize = 14)
plt.ylabel(r'$\xi^p$', fontsize = 14)
plt.ylim([-0.06,0.17])
plt.grid()
plt.gca().set_aspect('equal', adjustable='box')
plt.legend(loc=1)
plt.show()########