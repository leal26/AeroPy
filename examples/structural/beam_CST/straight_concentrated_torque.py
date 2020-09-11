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

g = CoordinateSystem.pCST(D=[0, 0, 0, 0, 0, 0], chord=[.5, .5], color=['b', 'g'],
                          N1=[1, 1], N2=[1, 1])
p = properties()
l = loads(torque=[1], torque_s=[.5])
s = np.linspace(0, 1, 101)

b = beam_chen(g, p, l, s, ignore_ends=True)
b.parameterized_solver(format_input, x0=g.D)
# b.g.calculate_x1(b.s)
b.calculate_resultants()
print(b.Rx, b.Ry)
print('Dm/ds', np.gradient(b.p.inertia*b.p.young*(b.g.rho-b.g_p.rho), b.s))
# b._residual([0.005714289063490557, 0.004285875208812112,
#              0.002857270578861466, -0.005714289063490557])
print('x', b.x)
print('y', b.y)
print('r', b.r)
print('R', b.R)
print('arc', b.g.darc)
print('rho', b.g.rho)
# Results from Abaqus
abaqus_data = pickle.load(open('neutral_line.p', 'rb'))

g0 = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord=1, color='k')
g0.calculate_x1(b.s)
g0.plot(label='Parent')

plt.scatter(abaqus[0, :], abaqus[1, :], c='.5', label='FEA', edgecolors='k', zorder=10, marker="^")
plt.plot(b.x, b.y, '.5', label='Child', lw=3, zorder=0)

plt.legend()
plt.show()
