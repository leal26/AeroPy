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
import math
import numpy as np
import warnings
from scipy.integrate import quad
from scipy.optimize import fsolve, minimize
from scipy import optimize
from scipy.optimize import differential_evolution

from airfoil_module import CST
from xfoil_module import output_reader

# Bersntein Polynomial
def K(r,n):
    K=math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
    return K

# Upper surface differential
def dxi_u(psi, Au, delta_xi):
    """Calculate upper derivate of xi for a given psi"""
    n = len(Au)-1
   
    diff = delta_xi/2.
    for i in range(n+1):
        diff += Au[i]*K(i,n)*(psi**i)*((1-psi)**(n-i))/(2*psi**0.5)*(-(3+2*n)*psi +2*i + 1)
    return  diff

# Lower surface differential
def dxi_l(psi, Al, delta_xi):
    """Calculate lower derivate of xi for a given psi"""
    n = len(Al)-1  
    diff = -delta_xi/2.
    for i in range(n+1):
        diff -= Al[i]*K(i,n)*psi**i*(1-psi)**(n-i)/(2*psi**0.5)*(-(3+2*n)*psi +2*i + 1)
    return diff

# Upper surface second differential
def ddxi_u(psi, Au, abs_output = False):
    """Calculate upper second derivate of xi for a given psi"""
    n = len(Au)-1
   
    diff = 0
    for i in range(n+1):
        diff -= Au[i]*K(i,n)*(psi**i)*((1-psi)**(n-i-1))/(4*psi**1.5)*((4*n**2 + 8*n +3)*psi**2+(-4*(2*i+1)*n - 4*i -2)*psi + 4*i**2 - 1)
    if abs_output:
        return abs(diff)
    else:
        return  diff

# Lower surface second differential
def ddxi_l(psi, Al, abs_output = False):
    """Calculate lower second derivate of xi for a given psi"""
    n = len(Al)-1  
    diff = 0
    for i in range(n+1):
        diff += Al[i]*K(i,n)*(psi**i)*((1-psi)**(n-i-1))/(4*psi**1.5)*((4*n**2 + 8*n +3)*psi**2+(-4*(2*i+1)*n-4*i-2)*psi +4*i**2 - 1)
    if abs_output:
        return abs(diff)
    else:
        return  diff
    
def calculate_c_baseline(c_L, Au_C, Au_L, deltaz):
    """Equations in the New_CST.pdf. Calculates the upper chord in order for
       the cruise and landing airfoils ot have the same length."""
    
    def integrand(psi, Au, delta_xi ):
        return np.sqrt(1 + dxi_u(psi, Au, delta_xi)**2)
    
    def f(c_C):
        """Function dependent of c_C and that outputs c_C."""
        y_C, err = quad(integrand, 0, 1, args=(Au_C, deltaz/c_C))
        y_L, err = quad(integrand, 0, 1, args=(Au_L, deltaz/c_L))
        return c_L*y_L/y_C
    c_C = optimize.fixed_point(f, [c_L])
    #In case the calculated chord is really close to the original, but the
    #algorithm was not able to make them equal
    if abs(c_L - c_C) < 1e-7:
        return c_L
    #The output is an array so it needs the extra [0]
    return c_C[0]

def calculate_psi_goal(psi_baseline, Au_baseline, Au_goal, deltaz,
                       c_baseline, c_goal):
    """Find the value for psi that has the same location w on the upper 
    surface of the goal as psi_baseline on the upper surface of the 
    baseline"""
    
    def integrand(psi_baseline, Au, deltaz, c ):
        return c*np.sqrt(1 + dxi_u(psi_baseline, Au, deltaz/c)**2)
    
    def equation(psi_goal, L_baseline, Au_goal, deltaz, c):
        y, err = quad(integrand, 0, psi_goal, args=(Au_goal, deltaz, c))
        return y - L_baseline
    
    L_baseline, err =  quad(integrand, 0, psi_baseline, args=(Au_baseline, deltaz, 
                                                         c_baseline))
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

def calculate_spar_direction(psi_baseline, Au_baseline, Au_goal, deltaz, c_goal):
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
        y_lower_goal = CST(psi_lower_goal*c_goal, c_goal, [deltaz/2., deltaz/2.], Au_goal, Al_goal)
        y_lower_goal = y_lower_goal['l']
        return psi_upper_goal + (s[0]/s[1])*(y_lower_goal - y_upper_goal)/c_goal
        
    # Calculate cruise chord
    c_baseline = calculate_c_baseline(c_goal, Au_baseline, Au_goal, deltaz)

    # Calculate upper psi at goal airfoil
    psi_upper_goal = calculate_psi_goal(psi_baseline, Au_baseline, Au_goal, 
                                        deltaz, c_baseline, c_goal)
    y_upper_goal = CST(psi_upper_goal*c_goal, c_goal, [deltaz/2., deltaz/2.], Au_goal, Al_goal)
    y_upper_goal = y_upper_goal['u']

    # Spar direction
    s = calculate_spar_direction(psi_baseline, Au_baseline, Au_goal, deltaz, c_goal)

    # Calculate lower psi and xi at goal airfoil
    #Because the iterative method can lead to warningdivision by zero after converging, we ignore
    #the warning
    np.seterr(divide='ignore', invalid='ignore') 
    psi_lower_goal = optimize.fixed_point(f, [psi_upper_goal]) #, args=(c_L, Au_C, Au_L, deltaz)
    x_lower_goal = psi_lower_goal*c_goal    
    y_lower_goal = CST(x_lower_goal, c_goal, [deltaz/2., deltaz/2.], Au_goal, Al_goal)
    y_lower_goal = y_lower_goal['l']

    return (y_upper_goal- y_lower_goal[0])/s[1]

def fitting_shape_coefficients(filename, bounds = 'Default', n = 5,
                               return_data = False, return_error = False,
                               optimize_deltaz = False):
    """Fit shape parameters to given data points
        Inputs:
        - filename: name of the file where the original data is
        - bounds: bounds for the shape parameters. If not defined,
                    Default values are used.
        - n: order of the Bernstein polynomial. If bounds is default
                this input will define the order of the polynomial.
                Otherwise the length of bounds (minus one) is taken into 
                consideration"""

    from hausdorff_distance import hausdorff_distance_2D

    def shape_difference(inputs, optimize_deltaz = False):

        if optimize_deltaz == True or optimize_deltaz == [True]:
            y_u = CST(upper['x'], 1, deltasz = inputs[-1]/2., Au = list(inputs[:n+1]))
            y_l = CST(lower['x'], 1, deltasz = inputs[-1]/2.,  Al = list(inputs[n+1:-1]))
        else:
            y_u = CST(upper['x'], 1, deltasz = deltaz/2., Au = list(inputs[:n+1]))
            y_l = CST(lower['x'], 1, deltasz = deltaz/2.,  Al = list(inputs[n+1:]))
        # Vector to be compared with
        a_u = {'x':upper['x'], 'y':y_u}
        a_l = {'x':lower['x'], 'y':y_l}
        
        b_u = upper
        b_l = lower
        return hausdorff_distance_2D(a_u, b_u) + hausdorff_distance_2D(a_l, b_l)

    # def shape_difference_upper(inputs, optimize_deltaz = False):
        # if optimize_deltaz == True:
            # y = CST(x, 1, deltasz = inputs[-1]/2., Au = list(inputs[:-1]))
        # else:
            # y = CST(x, 1, deltasz = inputs[-1]/2., Au = list(inputs))
        # # Vector to be compared with
        # b = {'x': x, 'y': y}
        # return hausdorff_distance_2D(a, b)

    # def shape_difference_lower(inputs, optimize_deltaz = False):
        # if optimize_deltaz == True:
            # y = CST(x, 1, deltasz = inputs[-1]/2.,  Al = list(inputs[:-1]))
        # else:
            # y = CST(x, 1, deltasz = deltaz/2.,  Al = list(inputs))
        # # Vector to be compared with
        # b = {'x': x, 'y': y}
        # return hausdorff_distance_2D(a, b)

    def separate_upper_lower(data):
        for i in range(len(data['x'])):
            if data['y'][i] < 0:
                break
        upper = {'x': data['x'][0:i],
                 'y': data['y'][0:i]}
        lower = {'x': data['x'][i:],
                 'y': data['y'][i:]}
        return upper, lower

    # Order of Bernstein polynomial
    if bounds != 'Default':
        n = len(bounds) - 1

    # Obtaining data
    data = output_reader(filename, separator = ', ', header = ['x', 'y'])

    # Rotating airfoil 
    x_TE = (data['x'][0] + data['x'][-1])/2.
    y_TE = (data['y'][0] + data['y'][-1])/2.

    theta_TE = math.atan(-y_TE/x_TE)

    # position trailing edge at the x-axis
    processed_data = {'x':[], 'y':[]}
    for i in range(len(data['x'])):
        x = data['x'][i]
        y = data['y'][i]
        c_theta = math.cos(theta_TE)
        s_theta = math.sin(theta_TE)
        x_rotated = c_theta*x - s_theta*y
        y_rotated = s_theta*x + c_theta*y
        processed_data['x'].append(x_rotated)
        processed_data['y'].append(y_rotated)
    data = processed_data

    # determine what is the leading edge and the rotation angle beta
    processed_data = {'x':[], 'y':[]}
    min_x_list = []
    min_y_list = []

    min_x = min(data['x'])
    min_index = data['x'].index(min_x)
    min_y = data['y'][min_index]

    chord = max(data['x']) - min(data['x'])
    beta = math.atan((y_TE - min_y)/(x_TE - min_x))
    
    for i in range(len(data['x'])):
        processed_data['x'].append((data['x'][i] - min_x)/chord)
        processed_data['y'].append(data['y'][i]/chord)    
    data = processed_data

    #==============================================================================
    # Optimizing shape
    #==============================================================================
    # Determining default bounds
    if bounds == 'Default':
        upper_bounds = [[0, 1.]]*(n+1)
        lower_bounds = [[0, 1]] +  [[-1., 1.]]*n

    if optimize_deltaz:
        bounds = upper_bounds + lower_bounds + [[0, 0.1]]
    else:
        bounds = upper_bounds + lower_bounds
        deltaz = (data['y'][0] - data['y'][-1])
    print bounds
    upper, lower = separate_upper_lower(data)
    # a = data
    # x = data['x']
    result = differential_evolution(shape_difference, bounds, 
                                            disp=True, popsize = 10, 
                                            args = [optimize_deltaz])
    print 'order %i upper done' % n
    # x = lower['x']
    # a = lower
    # result_lower = differential_evolution(shape_difference_lower, lower_bounds, 
                                            # disp=True, popsize = 10,
                                            # args = (optimize_deltaz))
    # print 'order %i lower done' % n
    if optimize_deltaz:
        Au = list(result.x[:n+1])
        Al = list(result.x[n+1:-1])
        deltaz = result.x[-1]
    else:
        Au = list(result.x[:n+1])
        Al = list(result.x[n+1:])
        
    # Return Al, Au, and others
    if return_data:
        return data, deltaz, Al, Au
    elif return_error:
        return result.fun, deltaz, Al, Au
    else:
        return deltaz, Al, Au

def shape_parameter_study(filename, n = 5):
    """Analyze the shape difference for different Bernstein order
       polynomials.
       - filename: name of dataset to compare with
       - n: Maximum Bernstein polynomial order """
    import pickle
    
    Data = {'error': [], 'Al': [], 'Au': [], 'order':[], 'deltaz':[]}
    for i in range(1,n+1):
        error, deltaz, Al, Au = fitting_shape_coefficients(filename, n = i, return_error = True, optimize_deltaz = True)
        Data['error'].append(error)
        Data['Al'].append(Al)
        Data['Au'].append(Au)
        Data['deltaz'].append(deltaz)
        Data['order'].append(i)
    
    file = open('shape_study.p', 'wb')
    pickle.dump(Data, file)
    return Data

def find_inflection_points(Au, Al):
    """Detect how many inflections points and where are they"""
    # Find solutions for several initial estimates
    x = np.linspace(0.000001,0.99999,100)
    psi_u_solutions = []
    psi_l_solutions = []
    
    # There will be solutions that do not converge, and will give
    # warnings. So just ignore them.
    # with warnings.catch_warnings():
        # warnings.simplefilter("ignore")
    for x_i in x:
        # Find solutions for upper and filter
        psi_i = minimize(ddxi_u, x_i, bounds =((0.0001,0.9999),), args=(Au, True))
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
        psi_i = minimize(ddxi_l, x_i, bounds =((0.0001,0.9999),), args=(Al, True))
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
    # ddxi = 0 is a necessary condition but not sufficient to be an inclination point
    # for such, a value right before and a value after need to have opposite signs
    m = len(psi_u_solutions)
    true_solutions_u = []
    psi_all = [0,] + list(psi_u_solutions) + [1,]

    if m != 0:
        for i in range(1, m+1):
            before_psi = (psi_all[i-1]+psi_all[i])/2.
            after_psi = (psi_all[i]+psi_all[i+1])/2.
            before_sign = np.sign(ddxi_u(before_psi,Au))
            after_sign = np.sign(ddxi_u(after_psi,Au))
            if after_sign + before_sign == 0:
                true_solutions_u.append(psi_all[i])

    m = len(psi_l_solutions)
    true_solutions_l = []
    psi_all = [0] + list(psi_l_solutions) + [1]
    if m != 0:
        for i in range(1, m+1):
            before_psi = (psi_all[i-1]+psi_all[i])/2.
            after_psi = (psi_all[i]+psi_all[i+1])/2.
            before_sign = np.sign(ddxi_u(before_psi,Au))
            after_sign = np.sign(ddxi_u(after_psi,Au))
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
    psi = np.linspace(0,1,1000)
    xi = CST(psi, 1., [delta_xi/2., delta_xi/2.], Au, Al)
    camber = (xi['u']+xi['l'])/2.
    return np.average(np.absolute(camber))
    
if __name__ == '__main__':
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt

    c_L = 0.36 # in meters
    deltaz = 0.002 #in meters
    
    Au_C = [0.4, 0.2]
    Au_L = [0.4, 0.2]
    
    Al_C = [0.4, 0.2]
    Al_L = [0.4, 0.2]

    c_C = calculate_c_baseline(c_L, Au_C, Au_L, deltaz)
    print "Solution:      ", c_C
    
    psi_i = 0.4
    print calculate_psi_goal(psi_i, Au_C, Au_L, deltaz, c_C, c_L)
    
    # Plot for several testing calculat_c_baseline
    x = np.linspace(0., 1., 11)
    
    print calculate_spar_distance(psi_i, Au_C, Au_L, Al_L, deltaz, c_L)
    c = []
    for x_i in x:
        Au_C[0] = x_i
        c_i = calculate_c_baseline(c_L, Au_C, Au_L, deltaz)
        c.append(c_i)
    plt.plot(x, c)
    plt.xlabel('$A_{u_0}^C$', fontsize = 20)
    plt.ylabel('$c^C$', fontsize = 20)
    plt.grid()
    plt.show()
    
    # Plot airfoils for different Au
    plt.figure()
    psi = np.linspace(0, 1, 500)
    i = 0
    for c_i in c:
        Au_C[0] = x[i]
        y = CST(psi, 1, [deltaz/2., deltaz/2.], Au_C, Al_C)
        x_plot = np.linspace(0, c_i, 500)
        plt.plot(x_plot, c_i*y['u'], label = '$A_{u_0}$ = %.1f' % x[i])
        y_psi = CST(psi_i, 1, [deltaz/2., deltaz/2.], Au_C, Al_C)
        i += 1
    plt.xlabel(r'$\psi^C$', fontsize = 20)
    plt.ylabel(r'$\xi^C$', fontsize = 20)
    plt.legend()
    plt.show()
    
    # Plot for several testing calculat_psi_goal
    plt.figure()
    x = np.linspace(0., 1.,11)
    psi_goal_list = []
    for x_i in x:
        Au_C[0] = x_i
        c_C = calculate_c_baseline(c_L, Au_C, Au_L, deltaz)
        psi_goal_i = calculate_psi_goal(psi_i, Au_C, Au_L, deltaz, c_C, c_L)
        psi_goal_list.append(psi_goal_i)
    plt.plot(x, psi_goal_list)
    plt.xlabel('$A_{u_0}^C$', fontsize = 20)
    plt.ylabel('$\psi_i^L$', fontsize = 20)
    plt.grid()
    plt.show()
    # Ploting psi_goal at the landing airfoil for different Au0 for cruise
    plt.figure()

    psi_plot = np.linspace(0, 1, 500)
    y = CST(psi_plot, 1, [deltaz/2., deltaz/2.], Au_L, Al_L)
    max_y = max(y['u'])
    plt.plot(psi_plot,y['u'])
   
    y = CST(psi_goal_list, 1, [deltaz/2., deltaz/2.], Au_L, Al_L)
    colors = iter(cm.rainbow(np.linspace(0, 1, len(psi_goal_list))))
    for i in range(len(psi_goal_list)):
        plt.scatter(psi_goal_list[i], y['u'][i], color=next(colors), label = '$A_{u_0}^C$ = %.1f' % x[i])
    plt.vlines(psi_i, 0, max_y, 'r')
    plt.xlabel('$\psi^L$', fontsize = 20)
    plt.ylabel(r'$\xi^L$', fontsize = 20)
    plt.legend()
    plt.show()

    # Plot cos(beta) several Au0
    plt.figure()
    x = np.linspace(0., 1.,11)
    cbeta_list = []
    for x_i in x:
        Au_C[0] = x_i
        cbeta_i = calculate_cbeta(psi_i, Au_C, deltaz/c_C)
        cbeta_list.append(cbeta_i)
        print Au_C, cbeta_i

    plt.plot(x, cbeta_list)
    plt.xlabel('$A_{u_0}^C$', fontsize = 20)
    plt.ylabel(r'cos($\beta$)', fontsize = 20)
    plt.grid()
    plt.show()
    
#    # Plot spar vector in landing and cruise configuration
    print 'Plot spar vector in landing and cruise configuration'
    x_landing = np.linspace(0, c_L, 500)
    plt.figure()
    colors = iter(cm.rainbow(np.linspace(0, 1, len(x))))
    for x_i in x:
        color_i = next(colors)
        Au_C[0] = x_i
        # Calculate cruise chord
        c_C = calculate_c_baseline(c_L, Au_C, Au_L, deltaz)
        x_cruise = np.linspace(0, c_C, 500)
        # Plot cruise airfoil
        y = CST(x_cruise, c_C, [deltaz/2., deltaz/2.], Au_C, Al_C)
        plt.plot(x_cruise, y['u'], c=color_i, label = '$A_{u_0}$ = %.2f' % x_i)
        plt.plot( x_cruise, y['l'], c=color_i)
        y = CST(psi_i*c_C, c_C, [deltaz/2., deltaz/2.], Au_C, Al_C)
        plt.plot([psi_i*c_C,psi_i*c_C], [y['l'], y['u']], c=color_i)
    plt.legend()
    plt.xlabel('$x^L$', fontsize = 20)
    plt.ylabel(r'$y^L$', fontsize = 20)
    plt.show()
    
    plt.figure()
    colors = iter(cm.rainbow(np.linspace(0, 1, len(x))))
    for x_i in x:
        Au_C[0] = x_i
        # Calculate cruise chord
        c_C = calculate_c_baseline(c_L, Au_C, Au_L, deltaz/c_L)
        # Calculate psi at landing
        psi_goal_i = calculate_psi_goal(psi_i, Au_C, Au_L, deltaz, c_C, c_L)
        # Calculate xi at landing
        temp = CST(psi_goal_i*c_L, c_L, [deltaz/2., deltaz/2.], Au_L, Al_L)
        y_goal_i = temp['u']
        # Plot landing airfoil
        y = CST(x_landing, c_L, [deltaz/2., deltaz/2.], Au_L, Al_L)
        plt.plot(x_landing, y['u'], 'b', x_landing, y['l'], 'b')
    
        #calculate spar direction
        s = calculate_spar_direction(psi_i, Au_C, Au_L, deltaz, c_L)
        #calculate spar length
        l = calculate_spar_distance(psi_i, Au_C, Au_L, Al_L, deltaz, c_L)
        print s, s[0]**2 + s[1]**2
        plt.scatter([psi_goal_i*c_L], [y_goal_i])
        plt.plot([psi_goal_i*c_L,psi_goal_i*c_L - l*s[0]],[y_goal_i, y_goal_i - l*s[1]], c = next(colors), label = '$A_{u_0}$ = %.2f' % x_i)
    plt.legend()
    plt.xlabel('$y^L$', fontsize = 20)
    plt.ylabel(r'$x^L$', fontsize = 20)
    plt.show()

#==============================================================================
#     #Plot to check if spars are the same if same airfoil
#==============================================================================
    plt.figure()
    Au_C = [0.2, 0.2]
    Au_L = [0.2, 0.2]
    
    Al_C = [0.2, 0.2]
    Al_L = [0.2, 0.2]
    
    c_C = calculate_c_baseline(c_L, Au_C, Au_L, deltaz)
    
    # Plot cruise airfoil
    x = np.linspace(0, c_C, 200)
    y = CST(x, c_C, deltasz= [deltaz/2., deltaz/2.],  Al= Al_C, Au =Au_C)
    plt.plot(x, y['u'], 'b', x, y['l'], 'b', label = 'landing')
    y = CST(psi_i*c_C, c_C, [deltaz/2., deltaz/2.], Au = Au_C, Al = Al_C)
    plt.plot([psi_i*c_C,psi_i*c_C], [y['l'], y['u']], 'b')
    
    print 'cruise spar y', y
    
    #Plot cruise spars
    x = np.linspace(0, c_L, 200)
    y = CST(x, c_L, deltasz= [deltaz/2., deltaz/2.],  Al= Al_L, Au =Au_L)
    plt.plot(x, y['u'], 'r', x, y['l'], 'r', label='cruise')
    
    #Find properties of landing spars
    # Calculate psi at landing
    psi_goal_i = calculate_psi_goal(psi_i, Au_C, Au_L, deltaz, c_C, c_L)
    x_goal_i = psi_goal_i*c_L
    # Calculate xi at landing
    temp = CST(x_goal_i, c_L, [deltaz/2., deltaz/2.], Au = Au_L, Al= Al_L)
    y_goal_i = temp['u']
    #calculate spar direction
    s = calculate_spar_direction(psi_i, Au_C, Au_L, deltaz, c_L)
    l = calculate_spar_distance(psi_i, Au_C, Au_L, Al_L, deltaz, c_L)
    plt.plot([x_goal_i, x_goal_i - l*s[0]],[y_goal_i, y_goal_i - l*s[1]], c = 'r')
    print 'landing spar y', [y_goal_i, y_goal_i - l*s[1]]
    plt.legend()
    plt.show()
#==============================================================================
#   Tests for curve fitting
#==============================================================================
    filename = 'sampled_airfoil_data.csv'
    # plt.figure()
    # # error, fitted_deltaz, fitted_Al, fitted_Au = fitting_shape_coefficients(filename, n=1,
                                                        # # optimize_deltaz = True, return_error = True)
    # # print error, fitted_Al, fitted_Au, fitted_deltaz
    # fitted_Au = [0.28112944407629581, 0.19487054006845383, 0.42397361498813563, 0.13538750382907982, 0.3399920533480057, 0.2118532593111192]
    # fitted_Al = [0.1897357530709628, -0.25258128225279725, 0.086871096306674597, -0.55958630302132484, 0.0064412971620611478, -0.24295645929089565]
    # fitted_deltaz = 0.0092804245707460483
    # data = output_reader(filename, separator = ', ', header = ['x', 'y'])
    # # Rotating airfoil 
    # x_TE = (data['x'][0] + data['x'][-1])/2.
    # y_TE = (data['y'][0] + data['y'][-1])/2.

    # theta_TE = math.atan(-y_TE/x_TE)

    # # position trailing edge at the x-axis
    # processed_data = {'x':[], 'y':[]}
    # for i in range(len(data['x'])):
        # x = data['x'][i]
        # y = data['y'][i]
        # c_theta = math.cos(theta_TE)
        # s_theta = math.sin(theta_TE)
        # x_rotated = c_theta*x - s_theta*y
        # y_rotated = s_theta*x + c_theta*y
        # processed_data['x'].append(x_rotated)
        # processed_data['y'].append(y_rotated)
    # data = processed_data

    # # plt.scatter(data['x'], data['y'])
    # x = np.linspace(0, 1, 200)
    # y = CST(x, 1, deltasz= [fitted_deltaz/2., fitted_deltaz/2.],  Al = fitted_Al, Au = fitted_Au)
    # plt.plot(x, y['u'], x, y['l'])
    # plt.scatter(data['x'], data['y'])
    # plt.show()
    # BREAK
#==============================================================================
#   Shape parameter study
#==============================================================================
    n = 8
    Data = shape_parameter_study(filename, n = n)
    plt.figure()
    x = np.linspace(2,2*n,n)
    plt.plot(x, Data['error'])
    plt.scatter(x, Data['error'])
    plt.grid()
    plt.xlabel('Number of shape functions')
    plt.ylabel('Haussdorf Distance (adimensional)')	
    plt.show()   