from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

import aeropy.xfoil_module as xf
from aeropy.CST.module_2D import *
from aeropy.aero_module import Reynolds
from aeropy.geometry.airfoil import CST, create_x

# Au = [0.23993240191629417, 0.34468227138908186, 0.18125405377549103, 
        # 0.35371349126072665, 0.2440815012119143, 0.25724974995738387]
# Al = [0.18889012559339036, -0.24686758992053115, 0.077569769493868401,
        # -0.547827192265256, -0.0047342206759065641, -0.23994805474814629]
Au= [0.172802, 0.167353, 0.130747, 0.172053, 0.112797, 0.168891]
Al = Au
# c_avian = 0.36                  #m
# deltaz = 0.0093943568219451313*c_avian
c_avian = 1.
deltaz = 0

airfoil = 'avian'
x = create_x(1., distribution = 'linear')
y = CST(x, 1., [deltaz/2., deltaz/2.], Au = Au, Al= Al)
# Create file for Xfoil to read coordinates
xf.create_input(x, y['u'], y['l'], airfoil, different_x_upper_lower = False)
print('Reynolds: ', Reynolds(10000, 30, c_avian))
Data = xf.find_coefficients(airfoil, 0., Reynolds=Reynolds(10000, 30, c_avian), iteration=100, NACA=False)
print(Data)