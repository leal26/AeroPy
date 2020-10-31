import matplotlib.pyplot as plt
import numpy as np
import math

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
n = 2
p = 3
i = (n+1)*p
g = CoordinateSystem.pCST(D=i*[0], chord=[.2, .4, .4], color=['b', 'g', 'r'],
                          N1=[1, 1, 1], N2=[1, 1, 1], continuity='C1', free_end=True,
                          root_fixed=True)
g.calculate_s(N=[21, 21, 21])
g.D = [2.4092723883736816e-07, 0.011416127507473277, 0.022835053534130918, 0.07336114552927285,
       0.0646572617461432, 0.055739427915895014, 0.036940013273059756, 0.02798017304135545]

g.calculate_x1(g.s)

for j in range(g.p):
    plt.figure()
    x = g.cst[j].x1_grid
    chord = g.cst[j].chord
    deltaT = chord*g.cst[j].zetaT
    deltaL = chord*g.cst[j].zetaL
    A = g.cst[j].D[:-1]
    for i in range(g.n+1):
        A_i = (g.n+1)*[0]
        A_i[i] = g.cst[j].D[i]
        g.cst[j].D[:-1] = A_i
        y = g.cst[j].x3(x, diff=diff)
        plt.plot(x, y, '--', label=i)
        g.cst[j].D[:-1] = A
    y = g.cst[j].x3(x, diff=diff)
    plt.plot(x, y, 'k')
    # plt.axis('equal')
    plt.grid()
    # plt.xlim([0, 1])
    plt.xlabel('$\psi$', fontsize=14)
    plt.ylabel(r'$\xi$', fontsize=14)
    plt.title('pCST %i' % (j+1))
    plt.legend()
plt.show()
