# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 20:59:22 2015

@author: Pedro
"""
from aero_module import *
from xfoil_module import *

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
    coefficients = find_coefficients(airfoil, alpha, Reynolds, iteration,
                                     NACA)
    alpha_L_0_root = find_alpha_L_0(airfoil, Reynolds, iteration, NACA)
    return LLT_calculator(alpha_L_0_root, coefficients['CD'], N, span, taper, chord_root,
                          alpha_root, velocity)


if __name__ == '__main__':
    print find_3D_coefficients(airfoil='naca0012', alpha=1.)
