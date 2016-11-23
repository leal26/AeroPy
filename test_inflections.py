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
c_avian = .36                  #m
deltaz = 0.0093943568219451313*c_avian

airfoil = 'avian'
x = create_x(c_avian, distribution = 'linear')
y = CST(x, c_avian, [deltaz/2., deltaz/2.], Au = Au, Al= Al)
# Create file for Xfoil to read coordinates
xf.create_input(x, y['u'], y['l'], airfoil, different_x_upper_lower = False)

Data = xf.find_coefficients(airfoil, 2., Reynolds=Reynolds(10000, 30, c_avian), iteration=100, NACA=False)
print Data
BREAK

psi_u_inflection, psi_l_inflection = find_inflection_points(Au, Al)
print 'upper: ', psi_u_inflection
print 'lower: ', psi_l_inflection
psi = np.linspace(0.001,0.999,100)
xi = CST(psi, 1, [deltaz/2., deltaz/2.], Au, Al)
plt.plot(psi, xi['u'], 'b', label = r'$\xi_u$')
plt.plot(psi, xi['l'],'b--', label = r'$\xi_l$')

xi_u_inflection = CST(psi_u_inflection, 1, [deltaz/2., deltaz/2.], Au, Al)
plt.scatter(psi_u_inflection, xi_u_inflection['u'])

xi_l_inflection = CST(psi_l_inflection, 1, [deltaz/2., deltaz/2.], Au, Al)
plt.scatter(psi_l_inflection, xi_l_inflection['l'])
plt.xlabel('$\psi$', fontsize = 20)
plt.ylabel(r'$\xi$', fontsize = 20)
plt.grid()
# plt.show()

# plt.figure()
plt.plot(psi,ddxi_u(psi,Au), 'r', label = r'$dd\xi_u$')
plt.plot(psi,ddxi_l(psi,Al), 'r--', label = r'$dd\xi_l$')
# plt.show()

# plt.figure()
plt.plot(psi,dxi_u(psi,Au, deltaz), 'g', label = r'$d\xi_u$')
plt.plot(psi,dxi_l(psi,Al,deltaz), 'g--', label = r'$d\xi_l$')
plt.ylim([-0.1,0.1])
plt.legend(loc='best')
plt.show()
