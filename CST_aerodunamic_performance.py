import numpy as np
import matplotlib.pyplot as plt

import xfoil_module as xf
from CST_module import *
from aero_module import Reynolds
from airfoil_module import CST, create_x

Au = [0.23993240191629417, 0.34468227138908186, 0.18125405377549103, 
        0.35371349126072665, 0.2440815012119143, 0.25724974995738387]
Al = [0.18889012559339036, -0.24686758992053115, 0.077569769493868401,
        -0.547827192265256, -0.0047342206759065641, -0.23994805474814629]
c_avian = 0.36                  #m
deltaz = 0.0093943568219451313*c_avian

airfoil = 'avian'
x = create_x(1., distribution = 'linear')
y = CST(x, 1., [deltaz/2., deltaz/2.], Au = Au, Al= Al)
# Create file for Xfoil to read coordinates
xf.create_input(x, y['u'], y['l'], airfoil, different_x_upper_lower = False)
print 'Reynolds: ', Reynolds(10000, 30, c_avian)
Data = xf.find_coefficients(airfoil, 2., Reynolds=Reynolds(10000, 30, c_avian), iteration=100, NACA=False)
print Data