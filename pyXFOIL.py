# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 20:59:22 2015

@author: Pedro
"""
from aero_tools import *
from xfoil_tools import *

def find_3D_coefficients(airfoil, alpha, Reynolds=0, iteration=10, NACA=True,
                         N=10, b=10., taper=1., chord_root=1, alpha_root=1.,
                         V=1.):
    coefficients = find_coefficients(airfoil, alpha, Reynolds, iteration,
                                     NACA)
    alpha_L_0_root = find_alpha_L_0(airfoil, Reynolds, iteration, NACA)
    return LLT_calculator(alpha_L_0_root, coefficients['CD'], N, b, taper, chord_root,
                          alpha_root, V)
                          

if __name__ == '__main__':
    print find_3D_coefficients(airfoil='naca0012', alpha=1.)