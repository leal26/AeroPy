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

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy import optimize

from airfoil_module import CST

def dxi_u(psi, Au0, Au1, delta_xi):
   """Calculate upper derivate of xi for a given psi"""
   return  (5*Au0*psi**2 - 3*(Au1 + 2*Au0)*psi + Au1 + Au0)/(2*np.sqrt(psi)) + \
           delta_xi/2.

def dxi_l(psi, Al0, Al1, delta_xi):
   """Calculate lower derivate of xi for a given psi"""
   return  - ((5*Al0*psi**2 - 3*(Al1 + 2*Al0)*psi + Al1 + Al0)/(2*np.sqrt(psi))
           + delta_xi/2.)
       
def calculate_c_baseline(c_L, Au_C, Au_L, deltaz):
    """Equations in the New_CST.pdf. Calculates the upper chord in order for
       the cruise and landing airfoils ot have the same length."""
    
    def integrand(psi, Au0, Au1, delta_xi ):
        return np.sqrt(1 + dxi_u(psi, Au0, Au1, delta_xi)**2)
    
    def f(c_C):
        """Function dependent of c_C and that outputs c_C."""
        y_C, err = quad(integrand, 0, 1, args=(Au_C[0], Au_C[1], deltaz/c_C))
        y_L, err = quad(integrand, 0, 1, args=(Au_L[0], Au_L[1], deltaz/c_L))
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
        return c*np.sqrt(1 + dxi_u(psi_baseline, Au[0], Au[1], deltaz/c)**2)
    
    def equation(psi_goal, L_baseline, Au_goal, deltaz, c):
        y, err = quad(integrand, 0, psi_goal, args=(Au_goal, deltaz, c))
        return y - L_baseline
    
    L_baseline, err =  quad(integrand, 0, psi_baseline, args=(Au_baseline, deltaz, 
                                                         c_baseline))

    y = fsolve(equation, psi_baseline, args=(L_baseline, Au_goal, deltaz,
                                                 c_goal))
    return y[0]
    
def calculate_cbeta(psi_i, Au0, Au1, delta_xi):
    """Calculate cosine for angles between vertical spars and outer mold
    line for cruise"""
    norm = np.sqrt(1+dxi_u(psi_i, Au0, Au1, delta_xi)**2)
    return dxi_u(psi_i, Au0, Au1, delta_xi)/norm

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
    
    cbeta = calculate_cbeta(psi_baseline, Au_baseline[0], Au_baseline[1], 
                            deltaz/c_baseline)
    sbeta = np.sqrt(1-cbeta**2)
    
    t[0] = 1
    t[1] = dxi_u(psi_goal, Au_goal[0], Au_goal[1], deltaz)
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
        return psi_upper_goal*c_goal + (s[0]/s[1])*(y_lower_goal - y_upper_goal)
        
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
    x_lower_goal = optimize.fixed_point(f, [psi_upper_goal]) #, args=(c_L, Au_C, Au_L, deltaz)
    y_lower_goal = CST(x_lower_goal, c_goal, [deltaz/2., deltaz/2.], Au_goal, Al_goal)
    y_lower_goal = y_lower_goal['l']

    return (y_upper_goal- y_lower_goal[0])/s[1]

if __name__ == '__main__':
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt

    c_L = 1. # in meters
    deltaz = 0.01 #in meters
    
    Au_C = [0.4, 0.3]
    Au_L = [0.2, 0.1]
    
    Al_C = [0.4, 0.3]
    Al_L = [0.2, 0.1]

    c_C = calculate_c_baseline(c_L, Au_C, Au_L, deltaz)
    print "Solution:      ", c_C
    
    psi_i = 0.2
    print calculate_psi_goal(psi_i, Au_C, Au_L, deltaz, c_C, c_L)
    
    
    # Plot for several testing calculat_c_baseline
    x = np.linspace(0., 1., 11)
    
    print calculate_spar_distance(psi_i, Au_C, Au_L, Al_L, deltaz, c_L)
#    c = []
#    for x_i in x:
#        Au_C[0] = x_i
#        c_i = calculate_c_baseline(c_L, Au_C, Au_L, deltaz)
#        c.append(c_i)
#    plt.plot(x, c)
#    plt.xlabel('$A_{u_0}^C$', fontsize = 20)
#    plt.ylabel('$c^C$', fontsize = 20)
#    plt.grid()
#    
#    # Plot airfoils for different Au
#    plt.figure()
#    psi = np.linspace(0, 1, 500)
#    i = 0
#    for c_i in c:
#        Au_C[0] = x[i]
#        y = CST(psi, 1, [deltaz/2., deltaz/2.], Au_C, Al_C)
#        x_plot = np.linspace(0, c_i, 500)
#        plt.plot(x_plot, c_i*y['u'], label = '$A_{u_0}$ = %.1f' % x[i])
#        y_psi = CST(psi_i, 1, [deltaz/2., deltaz/2.], Au_C, Al_C)
#        i += 1
#    plt.xlabel(r'$\psi^C$', fontsize = 20)
#    plt.ylabel(r'$\xi^C$', fontsize = 20)
#    plt.legend()
#    
#    # Plot for several testing calculat_psi_goal
#    plt.figure()
#    x = np.linspace(0., 1.,11)
#    psi_goal_list = []
#    for x_i in x:
#        Au_C[0] = x_i
#        c_C = calculate_c_baseline(c_L, Au_C, Au_L, deltaz)
#        psi_goal_i = calculate_psi_goal(psi_i, Au_C, Au_L, deltaz, c_C, c_L)
#        psi_goal_list.append(psi_goal_i)
#    plt.plot(x, psi_goal_list)
#    plt.xlabel('$A_{u_0}^C$', fontsize = 20)
#    plt.ylabel('$\psi_i^L$', fontsize = 20)
#    plt.grid()
#    
#    # Ploting psi_goal at the landing airfoil for different Au0 for cruise
#    plt.figure()
#
#    x_plot = np.linspace(0, c_L, 500)
#    y = CST(x_plot, 1, [deltaz/2., deltaz/2.], Au_L, Al_L)
#    max_y = max(y['u'])
#    plt.plot(x_plot,y['u'])
#    
#    y = CST(psi_goal_list, 1, [deltaz/2., deltaz/2.], Au_L, Al_L)
#    colors = iter(cm.rainbow(np.linspace(0, 1, len(psi_goal_list))))
#    for i in range(len(psi_goal_list)):
#        plt.scatter(psi_goal_list[i], y['u'][i], color=next(colors), label = '$A_{u_0}^C$ = %.1f' % x[i])
#    plt.vlines(psi_i, 0, max_y, 'r')
#    plt.xlabel('$\psi^L$', fontsize = 20)
#    plt.ylabel(r'$\xi^L$', fontsize = 20)
#    plt.legend()
#    
## Plot cos(beta) several Au0
#    plt.figure()
#    x = np.linspace(0., 1.,11)
#    cbeta_list = []
#    for x_i in x:
#        Au_C[0] = x_i
#        cbeta_i = calculate_cbeta(psi_i, Au_C[0], Au_C[1], deltaz/c_C)
#        cbeta_list.append(cbeta_i)
#    plt.plot(x, cbeta_list)
#    plt.xlabel('$A_{u_0}^C$', fontsize = 20)
#    plt.ylabel(r'cos($\beta$)', fontsize = 20)
#    plt.grid()
    
#    # Plot spar vector in landing and cruise configuration
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
    plt.figure()
#    
    colors = iter(cm.rainbow(np.linspace(0, 1, len(x))))
    for x_i in x:
        Au_C[0] = x_i
        # Calculate cruise chord
        c_C = calculate_c_baseline(c_L, Au_C, Au_L, deltaz)
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