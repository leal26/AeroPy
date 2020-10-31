import matplotlib.pyplot as plt
import numpy as np
import math
import copy

from aeropy.geometry.airfoil import create_x, CST
from aeropy.geometry.parametric import CoordinateSystem

diff = 'x11'
# Case study 1
# n = 2
# p = 2
# i = n*p+1
# g = CoordinateSystem.pCST(D=i*[0], chord=[.5, .5], color=['b', 'g'],
#                           N1=[1, 1], N2=[1, 1], continuity='C2', free_end=True,
#                           root_fixed=True)
# g.calculate_s(N=[11, 11])
# g.D = [4.078596014133988e-12, 0.0005361042099116627, 0.0010713304053242025, 0.0010708639962854962]

# Case study 2
g = CoordinateSystem.CST(D=[-0.27928065, - 0.29209946, - 0.30491828, - 0.3177371,
                            0.27928065], chord=0.4385, color='b', deltaz=0.27928065*0.4385, N1=1.0, N2=1.0)
g.calculate_s(N=21)
g.calculate_x1(g.s)
g_p = copy.copy(g)
# g.D = [-0.17142039268870685, -0.197078814282914, -
#        0.21528291934216998, -0.23359009242508835, 0.17142039268870685]

# g.calculate_x1(g.s)

plt.figure()
x = g.x1_grid
A = g.D[:-1]
for i in range(g.n+1):
    A_i = (g.n+1)*[0]
    A_i[i] = g.D[i]
    g.D[:-1] = A_i
    y = g.x3(x, diff=diff)
    plt.plot(x, y, '--', label=i)
    g.D[:-1] = A
y = g.x3(x, diff=diff)
plt.plot(x, y, 'k')
# plt.axis('equal')
plt.grid()
# plt.xlim([0, 1])
plt.xlabel('$\psi$', fontsize=14)
plt.ylabel(r'$\xi$', fontsize=14)
plt.legend()
plt.show()
