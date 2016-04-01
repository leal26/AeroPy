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

if __name__ == '__main__':
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
    Cm = calculate_flap_moment(Data['x'], Data['y'], alpha, x_hinge,
                               deflection, unit_deflection = 'rad')
    
    data = xf.find_pressure_coefficients('flapped', alpha, NACA = False)
    C_m = ar.calculate_moment_coefficient(data['x'], data['y'], data['Cp'], alpha)
    
    data_CM = xf.find_coefficients('flapped', alpha, NACA = False)
    print 'aeropy: overall: ', C_m
    print 'aeropy: flap: ', Cm
    
    print 'xfoil: overall: ', data_CM['CM']