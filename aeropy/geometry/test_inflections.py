from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

import aeropy.xfoil_module as xf
from aeropy.CST.module_2D import *
from aeropy.aero_module import Reynolds
from aeropy.airfoil_module import CST, create_x

Au = [0.23993240191629417, 0.34468227138908186, 0.18125405377549103, 
        0.35371349126072665, 0.2440815012119143, 0.25724974995738387]
Al = [0.18889012559339036, -0.24686758992053115, 0.077569769493868401,
        -0.547827192265256, -0.0047342206759065641, -0.23994805474814629]
c_avian = .36                  #m
deltaz = 0.0093943568219451313*c_avian

airfoil = 'avian'
x = create_x(c_avian, distribution = 'linear')
y = CST(x, c_avian, [deltaz/2., deltaz/2.], Au = Au, Al= Al)
# Create file for Xfoil to read coordinates
xf.create_input(x, y['u'], y['l'], airfoil, different_x_upper_lower = False)

Data = xf.find_coefficients(airfoil, 2., Reynolds=Reynolds(10000, 30, c_avian), iteration=100, NACA=False)
print(Data)


psi_u_inflection, psi_l_inflection = find_inflection_points(Au, Al)
print('upper: ', psi_u_inflection)
print('lower: ', psi_l_inflection)
psi = np.linspace(0.001,0.999,100)
xi = CST(psi, 1, [deltaz/2., deltaz/2.], Au, Al)
plt.plot(psi, xi['u'], 'b', label = 'Upper outer mold line')
plt.plot(psi, xi['l'],'b--', label = 'Lower outer mold line')

xi_u_inflection = CST(psi_u_inflection, 1, [deltaz/2., deltaz/2.], Au, Al)
plt.scatter(psi_u_inflection, xi_u_inflection['u'])

xi_l_inflection = CST(psi_l_inflection, 1, [deltaz/2., deltaz/2.], Au, Al)
plt.scatter(psi_l_inflection, xi_l_inflection['l'])
plt.xlabel('$\psi$', fontsize = 40)
plt.ylabel(r'$\xi$', fontsize = 40)
plt.grid()
# plt.show()

# plt.figure()
# plt.plot(psi,dxi_u(psi,Au, deltaz), 'g', label = r'$d\xi_u$')
# plt.plot(psi,dxi_l(psi,Al,deltaz), 'g--', label = r'$d\xi_l$')

plt.plot(psi,ddxi_u(psi,Au), 'r', label = 'Upper second derivative')
plt.plot(psi,ddxi_l(psi,Al), 'r--', label = 'Lower second derivative')

# Plot camber
camber = calculate_camber(psi, Au, Al, deltaz/c_avian)
plt.plot(psi,camber, 'k', label = 'camber')
psi_camber, xi_camber = calculate_max_camber(Au, Al, deltaz/c_avian)
print(psi_camber, xi_camber, type(psi_camber), type(xi_camber))
print('average camber: ', calculate_average_camber( Au, Al, deltaz/c_avian))
plt.scatter(psi_camber, xi_camber, label = 'Inflection points')
plt.ylim([-0.1,0.1])
plt.legend(loc='best')
plt.show()
