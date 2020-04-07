# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 14:41:37 2016

Suggested anity checks:
    - For same A coefficients, check if chords are equal
    - As A increases, the chord for cruise should shrink
    - For A coefficients, check if cos(beta)=0
    - For same A coefficients, check for calculate_psi_goal if psi_baseline is
      the same as psi_goal. Also try the initial guess from different points,
      to be sure it is good
@author: Pedro
"""
from __future__ import print_function
import math
import numpy as np
import warnings

from aeropy.geometry.airfoil import CST
from aeropy.xfoil_module import output_reader

try:
    from abaqus import *
    in_Abaqus = True
except(ModuleNotFoundError):
    in_Abaqus = False
    from scipy.integrate import quad
    from scipy.optimize import fsolve, minimize
    from scipy import optimize
    from scipy.optimize import differential_evolution


# Bersntein Polynomial
def K(r, n):
    K = math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
    return K

# Upper surface differential


def dxi_u(psi, Au, delta_xi, N1=0.5, N2=1):
    """Calculate upper derivate of xi for a given psi"""
    n = len(Au)-1

    xi_0 = CST(psi, 1, 0, Au, N1=N1, N2=N2)
    diff = xi_0*((1-n-N2))
    for i in range(n+1):
        # print N1-1., N2-1.
        # print psi**(N1-1.), (1-psi)**(N2-1.)
        # print Au[i]*K(i,n)*(psi**i)*((1-psi)**(n-i))*(i+N1-psi*(n+N1+N2))
        diff += (psi**(N1-1))*((1-psi)**(N2-1)) * \
            Au[i]*K(i, n)*(psi**i)*((1-psi)**(n-i))*(i+N1-psi*(n+N1+N2))
    return diff

# Lower surface differential


def dxi_l(psi, Al, delta_xi):
    """Calculate lower derivate of xi for a given psi"""
    n = len(Al)-1
    diff = -delta_xi/2.
    for i in range(n+1):
        diff -= Al[i] * K(i, n) * psi**i*(1-psi)**(n-i) / (2*psi**0.5) * \
            (-(3+2*n)*psi + 2*i + 1)
    return diff

# Upper surface second differential


def ddxi_u(psi, Au, abs_output=False):
    """Calculate upper second derivate of xi for a given psi"""
    n = len(Au)-1

    diff = 0
    for i in range(n+1):
        diff -= Au[i]*K(i, n)*(psi**i)*((1-psi)**(n-i-1))/(4*psi**1.5) * \
            ((4*n**2 + 8*n + 3)*psi**2+(-4*(2*i+1)*n - 4*i - 2)*psi +
             4*i**2 - 1)
    if abs_output:
        return abs(diff)
    else:
        return diff

# Lower surface second differential


def ddxi_l(psi, Al, abs_output=False):
    """Calculate lower second derivate of xi for a given psi"""
    n = len(Al)-1
    diff = 0
    for i in range(n+1):
        diff += Al[i]*K(i, n)*(psi**i)*((1-psi)**(n-i-1))/(4*psi**1.5) * \
            ((4*n**2 + 8*n + 3)*psi**2+(-4*(2*i+1)*n-4*i-2)*psi + 4*i**2 - 1)
    if abs_output:
        return abs(diff)
    else:
        return diff


def calculate_c_baseline(c_L, Au_C, Au_L, deltaz):
    """Equations in the New_CST.pdf. Calculates the upper chord in order for
       the cruise and landing airfoils ot have the same length."""

    def integrand(psi, Au, delta_xi):
        return np.sqrt(1 + dxi_u(psi, Au, delta_xi)**2)

    def f(c_C):
        """Function dependent of c_C and that outputs c_C."""
        y_C, err = quad(integrand, 0, 1, args=(Au_C, deltaz/c_C))
        y_L, err = quad(integrand, 0, 1, args=(Au_L, deltaz/c_L))
        return c_L*y_L/y_C
    c_C = optimize.fixed_point(f, [c_L])
    # In case the calculated chord is really close to the original, but the
    # algorithm was not able to make them equal
    if abs(c_L - c_C) < 1e-7:
        return c_L
    # The output is an array so it needs the extra [0]
    return c_C[0]


def calculate_psi_goal(psi_baseline, Au_baseline, Au_goal, deltaz,
                       c_baseline, c_goal):
    """Find the value for psi that has the same location w on the upper
    surface of the goal as psi_baseline on the upper surface of the
    baseline"""

    def integrand(psi_baseline, Au, deltaz, c):
        return c*np.sqrt(1 + dxi_u(psi_baseline, Au, deltaz/c)**2)

    def equation(psi_goal, L_baseline, Au_goal, deltaz, c):
        y, err = quad(integrand, 0, psi_goal, args=(Au_goal, deltaz, c))
        return y - L_baseline

    L_baseline, err = quad(integrand, 0, psi_baseline,
                           args=(Au_baseline, deltaz, c_baseline))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y = fsolve(equation, psi_baseline, args=(L_baseline, Au_goal, deltaz,
                                                 c_goal))
    return y[0]


def calculate_cbeta(psi_i, Au, delta_xi):
    """Calculate cosine for angles between vertical spars and outer mold
    line for cruise"""
    norm = np.sqrt(1+dxi_u(psi_i, Au, delta_xi)**2)
    return dxi_u(psi_i, Au, delta_xi)/norm


def calculate_spar_direction(psi_baseline, Au_baseline, Au_goal, deltaz,
                             c_goal):
    """Calculate the direction of the spar component based on a location
    at the upper surface for the cruise airfoil."""
    # Calculate cruise chord
    c_baseline = calculate_c_baseline(c_goal, Au_baseline, Au_goal, deltaz)
    # Calculate psi at goal arifoil
    psi_goal = calculate_psi_goal(psi_baseline, Au_baseline, Au_goal, deltaz,
                                  c_baseline, c_goal)
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


def calculate_spar_distance(psi_baseline, Au_baseline, Au_goal, Al_goal,
                            deltaz, c_goal):
    """Calculate spar distance (dimensional)"""

    def f(psi_lower_goal):
        y_lower_goal = CST(psi_lower_goal*c_goal, c_goal,
                           [deltaz/2., deltaz/2.], Au_goal, Al_goal)
        y_lower_goal = y_lower_goal['l']
        return psi_upper_goal + (s[0]/s[1])*(y_lower_goal -
                                             y_upper_goal)/c_goal

    # Calculate cruise chord
    c_baseline = calculate_c_baseline(c_goal, Au_baseline, Au_goal, deltaz)

    # Calculate upper psi at goal airfoil
    psi_upper_goal = calculate_psi_goal(psi_baseline, Au_baseline, Au_goal,
                                        deltaz, c_baseline, c_goal)
    y_upper_goal = CST(psi_upper_goal*c_goal, c_goal,
                       [deltaz/2., deltaz/2.], Au_goal, Al_goal)
    y_upper_goal = y_upper_goal['u']

    # Spar direction
    s = calculate_spar_direction(psi_baseline, Au_baseline, Au_goal, deltaz,
                                 c_goal)

    # Calculate lower psi and xi at goal airfoil
    # Because the iterative method can lead to warningdivision by zero after
    # converging, we ignore the warning
    np.seterr(divide='ignore', invalid='ignore')
    psi_lower_goal = optimize.fixed_point(f, [psi_upper_goal])
    x_lower_goal = psi_lower_goal*c_goal
    y_lower_goal = CST(x_lower_goal, c_goal, [deltaz/2., deltaz/2.],
                       Au_goal, Al_goal)
    y_lower_goal = y_lower_goal['l']

    return (y_upper_goal - y_lower_goal[0])/s[1]


def calculate_arc_length(psi_initial, psi_final, A_j, deltaz, c_j):
    """Calculate arc length from psi_initial to psi_final for
       shape coefficient A_j, trailing edge thickness deltaz, and
       chord c_j. Output is the dimensional length"""
    def integrand(psi_baseline, A_j, deltaz, c_j):
        return c_j*np.sqrt(1 + dxi_u(psi_baseline, A_j, deltaz/c_j)**2)

    L, err = quad(integrand, psi_initial, psi_final, args=(A_j, deltaz, c_j))
    return L


def find_inflection_points(Au, Al):
    """Detect how many inflections points and where are they"""
    # Find solutions for several initial estimates
    x = np.linspace(0.000001, 0.99999, 100)
    psi_u_solutions = []
    psi_l_solutions = []

    # There will be solutions that do not converge, and will give
    # warnings. So just ignore them.
    # with warnings.catch_warnings():
    # warnings.simplefilter("ignore")
    for x_i in x:
        # Find solutions for upper and filter
        psi_i = minimize(ddxi_u, x_i, bounds=((0.0001, 0.9999),),
                         args=(Au, True))
        psi_i = psi_i.x

        # Boolean to check if already in the list of solutions
        inside_upper = False
        for psi_j in psi_u_solutions:
            if psi_u_solutions == []:
                inside_upper = True
            elif abs(psi_j - psi_i[0]) < 1e-6:
                inside_upper = True
        # Boolean to check if actually a solution
        actual_solution = False

        if abs(ddxi_u(psi_i, Au))[0] < 1e-6:
            actual_solution = True
        if not inside_upper and actual_solution and (psi_i > 0 and psi_i < 1):
            psi_u_solutions.append(psi_i[0])

        # Find solutions for lower and filter
        psi_i = minimize(ddxi_l, x_i, bounds=((0.0001, 0.9999),),
                         args=(Al, True))
        psi_i = psi_i.x

        # Boolean to check if already in the list of solutions
        inside_lower = False
        for psi_j in psi_l_solutions:
            if psi_l_solutions == []:
                inside_lower = True
            elif abs(psi_j - psi_i[0]) < 1e-6:
                inside_lower = True
        # Boolean to check if actually a solution
        actual_solution = False

        if abs(ddxi_l(psi_i, Al))[0] < 1e-6:
            actual_solution = True
        if not inside_lower and actual_solution and (psi_i > 0 and psi_i < 1):
            psi_l_solutions.append(psi_i[0])
    # order lists
    psi_u_solutions = np.sort(psi_u_solutions)
    psi_l_solutions = np.sort(psi_l_solutions)
    # ddxi = 0 is a necessary condition but not sufficient to be an inclination
    # point for such, a value right before and a value after need to have
    # opposite signs
    m = len(psi_u_solutions)
    true_solutions_u = []
    psi_all = [0, ] + list(psi_u_solutions) + [1, ]

    if m != 0:
        for i in range(1, m+1):
            before_psi = (psi_all[i-1]+psi_all[i])/2.
            after_psi = (psi_all[i]+psi_all[i+1])/2.
            before_sign = np.sign(ddxi_u(before_psi, Au))
            after_sign = np.sign(ddxi_u(after_psi, Au))
            if after_sign + before_sign == 0:
                true_solutions_u.append(psi_all[i])

    m = len(psi_l_solutions)
    true_solutions_l = []
    psi_all = [0] + list(psi_l_solutions) + [1]
    if m != 0:
        for i in range(1, m+1):
            before_psi = (psi_all[i-1]+psi_all[i])/2.
            after_psi = (psi_all[i]+psi_all[i+1])/2.
            before_sign = np.sign(ddxi_u(before_psi, Au))
            after_sign = np.sign(ddxi_u(after_psi, Au))
            if after_sign + before_sign == 0:
                true_solutions_l.append(psi_all[i])
    return psi_u_solutions, psi_l_solutions


def calculate_camber(psi, Au, Al, delta_xi):
    xi = CST(psi, 1., [delta_xi/2., delta_xi/2.], Au, Al)
    return (xi['u']+xi['l'])/2.


def calculate_max_camber(Au, Al, delta_xi):
    """Calculate maximum camber and where it is. Returns (\psi, max_camber)"""
    def dcamber(psi, Au, Al, delta_xi):
        return 0.5*(dxi_u(psi, Au, delta_xi) + dxi_l(psi, Al, delta_xi))

    solution = fsolve(dcamber, 0.5, args=(Au, Al, delta_xi))

    # Outputs floats with psi and xi coordinates
    return solution[0], calculate_camber(solution, Au, Al, delta_xi)[0]


def calculate_average_camber(Au, Al, delta_xi):
    psi = np.linspace(0, 1, 1000)
    xi = CST(psi, 1., [delta_xi/2., delta_xi/2.], Au, Al)
    camber = (xi['u']+xi['l'])/2.
    return np.average(np.absolute(camber))
