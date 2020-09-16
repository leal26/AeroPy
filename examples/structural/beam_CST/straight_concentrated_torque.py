import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def format_input(input, g=None, g_p=None):
    # COnsidering BC for zero derivative at the root
    return list(input) + [-input[0]]
    # return input


abaqus = np.loadtxt('beam_torque.csv', delimiter=',')

g = CoordinateSystem.pCST(D=[0, 0, 0, 0, 0], chord=[.5, .5], color=['b', 'g'],
                          N1=[1, 1], N2=[1, 1], continuity='C1', free_end=True)

p = properties()
l = loads(torque=[1], torque_s=[.5])
g.calculate_s(N=[11, 11])
b = beam_chen(g, p, l, s=None, ignore_ends=True)
b.parameterized_solver(format_input, x0=g.D[:-1])
# [-0.00428542530710673, -0.004284305844033915, -0.00428524698261122, -
#                                          3.5282148257309257e-06, -3.054801179069392e-07, -2.9282143007525176e-06, 0.0, 0.00428542530710673]
# b.g.calculate_x1(b.s)
# b.calculate_resultants()
# print(b.Rx, b.Ry)
plt.figure()
plt.plot(b.s, b.p.inertia*b.p.young*(b.g.rho-b.g_p.rho), c='r')
print('s', b.s)
print('M', b.M)
plt.plot(b.s, b.M, c='b')
# b._residual([0.005714289063490557, 0.004285875208812112,
#              0.002857270578861466, -0.005714289063490557])
# print('x', b.x)
# print('y', b.y)
# print('r', b.r)
# print('R', b.R)
# print('arc', b.g.darc)
# print('rho', b.g.rho)
print('states')
print(b.g.n)
print(b.g.cst[0].zetaT, b.g.cst[1].zetaT)
print(b.g.cst[0].zetaL, b.g.cst[1].zetaL)
print(b.g.cst[0].chord, b.g.cst[1].chord)
print(b.g.cst[0].D, b.g.cst[1].D)
print()
# Results from Abaqus
abaqus_data = pickle.load(open('neutral_line.p', 'rb'))

g0 = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord=1, color='k')
g0.calculate_x1(b.s)
plt.figure()
g0.plot(label='Parent')

plt.scatter(abaqus[0, :], abaqus[1, :], c='.5', label='FEA', edgecolors='k', zorder=10, marker="^")
plt.plot(b.x, b.y, '.5', label='Child', lw=3, zorder=0)

plt.legend()
plt.show()
