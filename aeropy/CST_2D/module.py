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

    if N1 == 0.5 and N2 == 1:
        diff = delta_xi/2.
    else:
        diff = delta_xi
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


def fitting_shape_coefficients(filename, bounds='Default', n=5,
                               return_data=False, return_error=False,
                               optimize_deltaz=False, solver='gradient',
                               deltaz=None, objective='hausdorf',
                               surface='both'):
    """Fit shape parameters to given data points
        Inputs:
        - filename: name of the file where the original data is
        - bounds: bounds for the shape parameters. If not defined,
                    Default values are used.
        - n: order of the Bernstein polynomial. If bounds is default
                this input will define the order of the polynomial.
                Otherwise the length of bounds (minus one) is taken into
                consideration"""

    from optimization_tools.hausdorff_distance import hausdorff_distance_2D

    def shape_difference(inputs, optimize_deltaz=False, surface=surface):
        # Define deltaz
        if optimize_deltaz is True or optimize_deltaz == [True]:
            deltasz = inputs[-1]/2.
        else:
            deltasz = deltaz/2.

        # Calculate upper and lower surface
        if surface == 'both':
            y_u = CST(upper['x'], 1, deltasz=inputs[-1]/2.,
                      Au=list(inputs[:n+1]))
            y_l = CST(lower['x'], 1, deltasz=inputs[-1]/2.,
                      Al=list(inputs[n+1:-1]))
        elif surface == 'upper':
            y_u = CST(upper['x'], 1, deltasz=inputs[-1]/2.,
                      Au=list(inputs[:n+1]))
        elif surface == 'lower':
            y_l = CST(lower['x'], 1, deltasz=inputs[-1]/2.,
                      Al=list(inputs[:n+1]))

        # Vector to be compared with
        error = 0
        if surface == 'upper' or surface == 'both':
            a_u = {'x': upper['x'], 'y': y_u}
            if objective == 'hausdorf':
                error += hausdorff_distance_2D(a_u, upper)
            elif objective == 'squared_mean':
                error += np.mean((np.array(a_u['x'])-np.array(upper['x']))**2 +
                                 (np.array(a_u['y'])-np.array(upper['y']))**2)

        if surface == 'lower' or surface == 'both':
            a_l = {'x': lower['x'], 'y': y_l}
            if objective == 'hausdorf':
                error += hausdorff_distance_2D(a_l, lower)
            elif objective == 'squared_mean':
                error += np.mean((np.array(a_l['x'])-np.array(lower['x']))**2 +
                                 (np.array(a_l['y'])-np.array(lower['y']))**2)

        # plt.figure()
        # plt.scatter(a_u['x'], a_u['y'], c='k')
        # plt.scatter(a_l['x'], a_l['y'], c='b')
        # plt.scatter(upper['x'], upper['y'], c='r')
        # plt.scatter(lower['x'], lower['y'], c='g')
        # plt.show()
        return error

    def separate_upper_lower(data):
        for key in data:
            data[key] = np.array(data[key])

        index = np.where(data['y'] > 0)
        upper = {'x': data['x'][index],
                 'y': data['y'][index]}
        index = np.where(data['y'] <= 0)
        lower = {'x': data['x'][index],
                 'y': data['y'][index]}
        # x = data['x']
        # y = data['y']
        # for i in range(len(x)):
        #     if data['y'][i] < 0:
        #         break
        # upper = {'x': x[0:i],
        #          'y': y[0:i]}
        # lower = {'x': x[i:],
        #          'y': y[i:]}
        return upper, lower

    def _determine_bounds_x0(n, optimize_deltaz, bounds):
        if bounds == 'Default':
            upper_bounds = [[0, 1]] + [[-1., 1.]]*n
            lower_bounds = [[0, 1]] + [[-1., 1.]]*n

        if optimize_deltaz:
            if surface == 'both':
                bounds = upper_bounds + lower_bounds + [[0, 0.1]]
                x0 = (n+1)*[0., ] + (n+1)*[0., ] + [0.]
            elif surface == 'upper':
                bounds = upper_bounds + [[0, 0.1]]
                x0 = (n+1)*[0., ] + [0.]
            elif surface == 'lower':
                bounds = lower_bounds + [[0, 0.1]]
                x0 = (n+1)*[0., ] + [0.]
        else:
            bounds = upper_bounds + lower_bounds
            x0 = (n+1)*[0., ] + (n+1)*[0., ]
        return x0, bounds

    # Order of Bernstein polynomial
    if bounds != 'Default':
        n = len(bounds) - 1

    # Obtaining data
    if filename[-2:] == '.p':
        import pickle
        data = pickle.load(open(filename, "rb"), encoding='latin1')
        data = data['wing'][list(data['wing'].keys())[3]]
        x, y, z = data.T
    else:
        data = output_reader(filename, separator='\t', header=['x', 'z'])
        x = data['x']
        z = data['z']

    # Rotating airfoil
    x_TE = (x[0] + x[-1])/2.
    y_TE = (z[0] + z[-1])/2.

    theta_TE = math.atan(-y_TE/x_TE)

    # position trailing edge at the x-axis
    processed_data = {'x': [], 'y': []}
    for i in range(len(x)):
        x_i = x[i]
        z_i = z[i]
        c_theta = math.cos(theta_TE)
        s_theta = math.sin(theta_TE)
        x_rotated = c_theta*x_i - s_theta*z_i
        z_rotated = s_theta*x_i + c_theta*z_i
        processed_data['x'].append(x_rotated)
        processed_data['y'].append(z_rotated)
    data = processed_data

    # determine what is the leading edge and the rotation angle beta
    processed_data = {'x': [], 'y': []}

    min_x = min(x)
    # min_index = data['x'].index(min_x)
    # min_y = data['y'][min_index]

    chord = max(x) - min(x)
    # beta = math.atan((y_TE - min_y)/(x_TE - min_x))

    for i in range(len(x)):
        processed_data['x'].append((x[i] - min_x)/chord)
        processed_data['y'].append(z[i]/chord)
    data = processed_data

    # Determining default bounds
    x0, bounds = _determine_bounds_x0(n, optimize_deltaz, bounds)

    if not optimize_deltaz and deltaz is None:
        deltaz = (data['y'][0] - data['y'][-1])

    if surface == 'both':
        upper, lower = separate_upper_lower(data)
    elif surface == 'upper':
        upper = data
    elif surface == 'lower':
        lower = data

    # Calculate original error
    error0 = shape_difference(x0, optimize_deltaz=optimize_deltaz,
                              surface=surface)

    def f(x):
        return shape_difference(x, optimize_deltaz=optimize_deltaz,
                                surface=surface)/error0
    # Optimize
    if solver == 'differential_evolution':

        result = differential_evolution(f, bounds,
                                        disp=True, popsize=10)
        x = result.x
        f = result.fun
    elif solver == 'gradient':

        solution = minimize(f, x0, bounds=bounds,
                            options={'maxfun': 30000, 'eps': 1e-02})
        x = solution['x']
        f = solution['fun']
    print('order %i  done' % n)

    # Unpackage data
    if surface == 'both' or surface == 'upper':
        Au = list(x[:n+1])
    if surface == 'both':
        if optimize_deltaz:
            Al = list(x[n+1:-1])
            deltaz = x[-1]
        else:
            Al = list(x[n+1:])
    elif surface == 'lower':
        Al = list(x[:n+1])

    # Return Al, Au, and others
    to_return = []
    if return_data:
        to_return.append(data)
    if return_error:
        to_return.append(f)
    to_return.append(deltaz)
    if surface == 'lower' or surface == 'both':
        to_return.append(Al)
    elif surface == 'upper' or surface == 'both':
        to_return.append(Au)
    print(to_return)
    return to_return


def shape_parameter_study(filename, n=5, solver='gradient', deltaz=None,
                          objective='hausdorf', surface='both'):
    """Analyze the shape difference for different Bernstein order
       polynomials.
       - filename: name of dataset to compare with
       - n: Maximum Bernstein polynomial order """
    import pickle

    if deltaz is None:
        optimize_deltaz = True
    else:
        optimize_deltaz = False
    Data = {'error': [], 'Al': [], 'Au': [], 'order': [], 'deltaz': []}
    for i in range(1, n+1):
        if surface == 'both':
            error, deltaz, Al, Au = fitting_shape_coefficients(
                filename, n=i, return_error=True,
                optimize_deltaz=optimize_deltaz, solver=solver, deltaz=deltaz,
                objective=objective, surface=surface)
        elif surface == 'lower':
            error, deltaz, Al = fitting_shape_coefficients(
                filename, n=i, return_error=True,
                optimize_deltaz=optimize_deltaz, solver=solver, deltaz=deltaz,
                objective=objective, surface=surface)
        elif surface == 'upper':
            error, deltaz, Au = fitting_shape_coefficients(
                filename, n=i, return_error=True,
                optimize_deltaz=optimize_deltaz, solver=solver, deltaz=deltaz,
                objective=objective, surface=surface)
        print(error)
        Data['error'].append(error)
        if surface == 'both' or surface == 'lower':
            Data['Al'].append(Al)
        if surface == 'both' or surface == 'upper':
            Data['Au'].append(Au)
        Data['deltaz'].append(deltaz)
        Data['order'].append(i)

    file = open('shape_study.p', 'wb')
    pickle.dump(Data, file)
    return Data


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


if __name__ == '__main__':
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt

# ==============================================================================
#   Tests for curve fitting
# ==============================================================================
    # plt.figure()
    filename = '../filehandling/examples/SCFCoordinates.txt'
    surface = 'lower'
    output = fitting_shape_coefficients(filename, n=8, optimize_deltaz=False,
                                        return_error=True, deltaz=0,
                                        solver='differential_evolution',
                                        objective='squared_mean',
                                        surface=surface)
    if surface == 'both':
        error, fitted_deltaz, fitted_Al, fitted_Au = output
        print(error)
        print(fitted_Al)
        print(fitted_Au)
    elif surface == 'lower':
        error, fitted_deltaz, fitted_Al = output

    data = output_reader(filename, separator='\t', header=['x', 'z'])
    plt.scatter(data['x'], data['z'], c='r')
    x = np.linspace(0, 1, 100)
    # y_u = CST(x, 1, deltasz=0, Au=fitted_Au)
    y_l = CST(x, 1, deltasz=0, Al=fitted_Al)
    # plt.plot(x, y_u, 'b')
    plt.plot(x, y_l, 'b')
    plt.show()

# ==============================================================================
#   Shape parameter study
# ==============================================================================
    n = 8
    Data = shape_parameter_study(filename, n=n, solver='gradient', deltaz=0,
                                 objective='squared_mean', surface='lower')
    plt.figure()
    x = np.linspace(2, 2*n, n)
    plt.plot(x, Data['error'])
    plt.scatter(x, Data['error'])
    plt.grid()
    plt.xlabel('Number of shape functions')
    plt.ylabel('Haussdorf Distance (adimensional)')
    plt.show()
