# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 20:59:22 2015

@author: Pedro
"""
import aero_module as ar
import xfoil_module as xf
import airfoil_module as af

def find_3D_coefficients(airfoil, alpha, Reynolds=0, iteration=10, NACA=True,
                         N=10, span=10., taper=1., chord_root=1, alpha_root=1.,
                         velocity=1.):
    """ Calculate the 3D distribution using the Lifting Line Theory.
    
    :param airfoil: if NACA is false, airfoil is the name of the plain
           filewhere the airfoil geometry is stored (variable airfoil).
           If NACA is True, airfoil is the naca series of the airfoil
           (i.e.: naca2244). By default NACA is False.

    :param Reynolds: Reynolds number in case the simulation is for a
          viscous flow. In case not informed, the code will assume
          inviscid. (Use the aero_module function to calculate reynolds)
          
    :param alpha: list/array/float/int of angles of attack.

    :param iteration: changes how many times XFOIL will try to make the
          results converge. Specialy important for viscous flows

    :param NACA: Boolean variable that defines if the code imports an
          airfoil from a file or generates a NACA airfoil.
    
    :param N: number of cross sections on the wing
    
    :param span: span in meters
    
    :param taper: unidimendional taper (This options is still not 100%
            operational)
    
    :param chord_root: value of the chord at the the root
    
    :param alpha_root: angle of attack of the chord at the root (degrees)

    :param velocity: velocity in m/s

"""
    coefficients = xf.find_coefficients(airfoil, alpha, Reynolds, iteration,
                                     NACA)
    alpha_L_0_root = xf.find_alpha_L_0(airfoil, Reynolds, iteration, NACA)
    return ar.LLT_calculator(alpha_L_0_root, coefficients['CD'], N, span, taper, chord_root,
                          alpha_root, velocity)

def calculate_flap_moment(x, y, alpha, x_hinge, deflection,
                          unit_deflection = 'rad'):    
    """For a given airfoil with coordinates x and y at angle of attack
    alpha (degrees), calculate the moment coefficient around the joint at x_hinge
    and deflection in radians (unit_deflection = 'rad') or degrees 
    (unit_deflection = 'deg')"""

    # If x and y are not dictionaries with keys upper and lower, make them
    # be so
    if type(x) == list:
        x, y = af.separate_upper_lower(x, y)
    #Because parts of the program use keys 'u' and 'l', and other parts use
    #'upper' and 'lower'
    if 'u' in x.keys():
        upper = {'x': x['u'], 'y': y['u']}
        lower = {'x': x['l'], 'y': y['l']}    
    elif 'upper' in x.keys(): 
        upper = {'x': x['upper'], 'y': y['upper']}
        lower = {'x': x['lower'], 'y': y['lower']}
    
    hinge = af.find_hinge(x_hinge, upper, lower)
       
    if deflection > 0:
        upper_static, upper_flap = af.find_flap(upper, hinge)
        lower_static, lower_flap = af.find_flap(lower, hinge,
                                                extra_points = 'lower')
    elif deflection < 0:
        upper_static, upper_flap = af.find_flap(upper, hinge,
                                                extra_points = 'upper')
        lower_static, lower_flap = af.find_flap(lower, hinge)
    else:
       upper_static, upper_flap = af.find_flap(upper, hinge, 
                                               extra_points = None)
       lower_static, lower_flap = af.find_flap(lower, hinge,
                                               extra_points = None)

    upper_rotated, lower_rotated = af.rotate(upper_flap, lower_flap,
                                             hinge, deflection,
                                             unit_theta = unit_deflection)
     
    flapped_airfoil, i_separator = af.clean(upper_static, upper_rotated, lower_static, 
                            lower_rotated, hinge, deflection, N = None, 
                            return_flap_i = True)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Third step: get new values of pressure coefficient
    xf.create_input(x = flapped_airfoil['x'], y_u = flapped_airfoil['y'],
                    filename = 'flapped', different_x_upper_lower = True)

    Data = xf.find_pressure_coefficients('flapped', alpha, NACA = False)


    x, y, Cp = af.separate_upper_lower(x = Data['x'], y = Data['y'], 
                                        Cp = Data['Cp'], i_separator = i_separator)
    
    #At the hinges, the reaction moment has the opposite sign of the
    #actuated torque
    Cm = - ar.calculate_moment_coefficient(x, y, Cp, alpha = alpha, c = 1., 
                                         x_ref = x_hinge, y_ref = 0., 
                                         flap = True)
    return Cm

def calculate_flap_coefficients(x, y, alpha, x_hinge, deflection, Reynolds = 0,):    
    """For a given airfoil with coordinates x and y at angle of attack
    alpha, calculate the moment coefficient around the joint at x_hinge
    and deflection """

    def separate_upper_lower(x, y, Cp = None, i_separator=None):
        """Return dictionaries with upper and lower keys with respective
        coordiantes. It is assumed the leading edge is frontmost point at
        alpha=0"""
        #TODO: when using list generated by xfoil, there are two points for
        #the leading edge
        def separate(variable_list, i_separator):
            if type(i_separator) == int:
                variable_dictionary = {'upper': variable_list[0:i_separator+1],
                                       'lower': variable_list[i_separator+1:]}
            elif type(i_separator) == list:
                i_upper = i_separator[0]
                i_lower = i_separator[1]
                
                variable_dictionary = {'upper': variable_list[0:i_upper],
                                       'lower': variable_list[i_lower:]}
            return variable_dictionary
        #If i is not defined, separate upper and lower surface from the
        # leading edge
        if i_separator == None:
            i_separator = x.index(min(x))

        if Cp == None:
            x = separate(x, i_separator)
            y = separate(y, i_separator)
            return x, y
        else:
            x = separate(x, i_separator)
            y = separate(y, i_separator)
            Cp = separate(Cp, i_separator)
            return x, y, Cp
    # If x and y are not dictionaries with keys upper and lower, make them
    # be so
    if type(x) == list:
        x, y = separate_upper_lower(x, y)
    
    upper = {'x': x['upper'], 'y': y['upper']}
    lower = {'x': x['lower'], 'y': y['lower']}
    
    #Determining hinge
    hinge = af.find_hinge(x_hinge, upper, lower)
    
    upper_static, upper_flap = af.find_flap(upper, hinge)
    lower_static, lower_flap = af.find_flap(lower, hinge, lower = True)
    
    upper_rotated, lower_rotated = af.rotate(upper_flap, lower_flap, hinge, deflection)
    
    flapped_airfoil, i_separator = af.clean(upper_static, upper_rotated, lower_static, 
                            lower_rotated, hinge, deflection, N = 5, 
                            return_flap_i = True)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Third step: get new values of pressure coefficient
    xf.create_input(x = flapped_airfoil['x'], y_u = flapped_airfoil['y'],
                    filename = 'flapped', different_x_upper_lower = True)
    Data = xf.find_coefficients('flapped', alpha, NACA = False,
                                Reynolds = Reynolds, iteration=500)
    
    return Data


if __name__ == '__main__':

    import numpy as np
    import matplotlib.pyplot as plt
    import math

#    print find_3D_coefficients(airfoil='naca0012', alpha=1.)
    alpha = 0.
    x_hinge = 0.25/0.6175
    deflection = -math.pi/2. #0.17453292519943295 #0.0010573527055
    
    # generate original airfoil
    airfoil = "naca0012"
    xf.call(airfoil, output='Coordinates')
    filename = xf.file_name(airfoil, output='Coordinates')
    Data = xf.output_reader(filename, output='Coordinates', header = ['x','y'])

    Cm = calculate_flap_moment(Data['x'], Data['y'], alpha, x_hinge, deflection)

    
    V = 10
    altitude = 10000 #feet
    
    Reynolds = ar.Reynolds(altitude, V, 1.0)
    
    deflection_list = list(np.linspace(5,30,4))
    alpha_list = list(np.linspace(0,15,20))

    # Calculate coefficients for without flap
    CL_list = []
    CD_list = []
    ratio_list = []
    CM_list = []
        
    for alpha_i in alpha_list:
        Data_0 = xf.find_coefficients('naca0012', alpha_i, Reynolds = Reynolds,
                                      iteration = 200)
        CL_list.append(Data_0['CL'])
        CD_list.append(Data_0['CD'])
        ratio_list.append(Data_0['CL']/Data_0['CD'])
        CM_list.append(Data_0['CM'])    
    All_data = {0:{r'$c_m$': CM_list, r'$c_l$':CL_list,
                     r'$c_d$': CD_list,
                     r'$c_l/c_d$': ratio_list}}#:Data_0['CL']/Data_0['CD']}}
                     
    # Calculate foeccifient when using flap
    for deflection_i in deflection_list:
        CL_list = []
        CD_list = []
        ratio_list = []
        CM_list = []
        for alpha_i in alpha_list: 
            flap_data = calculate_flap_coefficients(Data['x'], Data['y'], alpha_i, x_hinge,
                                               deflection_i, Reynolds = Reynolds)
            CL_list.append(flap_data['CL'])
            CD_list.append(flap_data['CD'])
            ratio_list.append(flap_data['CL']/flap_data['CD'])
            CM_list.append(flap_data['CM'])
        All_data[deflection_i] = {r'$c_m$' : CM_list, r'$c_l$': CL_list,
                                  r'$c_d$' : CD_list, r'$c_l/c_d$': ratio_list}
    for key in [r'$c_m$', r'$c_l$', r'$c_d$', r'$c_l/c_d$']:
        plt.figure()
        for deflection_i in [0] + deflection_list:
            plt.plot(alpha_list, All_data[deflection_i][key], label = r'$\theta$ = %.0f' % deflection_i)
        plt.legend(loc = "best")
        plt.xlabel(r'$\alpha$', fontsize = 22)
        plt.ylabel(key, fontsize = 22)
        plt.grid()
 
    plt.figure()
    for deflection_i in [0] + deflection_list:
        plt.plot(All_data[deflection_i][r'$c_d$'], 
                 All_data[deflection_i][r'$c_l$'], 
                 label = r'$\theta$ = %.0f' % deflection_i)
    plt.legend(loc = "best")
    plt.xlabel(r'$c_d$', fontsize = 22)
    plt.ylabel(r'$c_l$', fontsize = 22)
    plt.grid()               
