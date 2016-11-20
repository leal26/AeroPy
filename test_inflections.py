import numpy as np
import matplotlib.pyplot as plt

from CST_module import *
from airfoil_module import CST

Au = [0.2, -0.1, 0.2, -0.1, 0.2]
Al = [0.2, -0.1, 0.2, -0.1, 0.2]
deltaz = 0.001

psi_u_inflection, psi_l_inflection, lengths = find_inflection_points(Au, Al)
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
