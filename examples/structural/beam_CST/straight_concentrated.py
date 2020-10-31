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


g = CoordinateSystem.CST(D=[0, 0, 0, 0], chord=1, color='b', N1=1, N2=1)
p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[1])
s = np.linspace(0, 1, 10)

b = beam_chen(g, p, l, s, ignore_ends=True)
b.parameterized_solver(format_input, x0=g.D[:-1])

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

# Euler-Bernoulle
EB_solution = l.concentrated_load[0][-1]/(6*p.young*p.inertia) * \
    np.array([0, 0, 3, -1])
print('EB', EB_solution)
curve_EB = CoordinateSystem.polynomial(D=EB_solution, chord=1, color='0.5')
curve_EB.calculate_x1(b.s)
eulerBernoulle = curve_EB.r(b.s)

g0 = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord=1, color='k')
g0.calculate_x1(b.s)
g0.plot(label='Parent')
[x, y] = eulerBernoulle.T
print('EB x:', x)
print('EB y:', y)
plt.scatter(x, y, c='.5', label='Euler-Bernoulle: %.3f N' % -
            l.concentrated_load[0][1], linestyle='-', zorder=1, edgecolors='k')
plt.scatter(abaqus_data['coord'][0:401:40, 0],
            abaqus_data['coord'][0:401:40, 1], c='.5',
            label='FEA: %.3f N' % -l.concentrated_load[0][1],
            edgecolors='k', zorder=2, marker="^")
plt.plot(b.x, b.y, '.5', label='Child: %.3f N' % -l.concentrated_load[0][1],
         lw=3, zorder=0)

plt.legend()

plt.figure()
plt.plot(b.g.x1_grid, b.g.rho)

plt.figure()
plt.plot(b.g.x1_grid[1:], b.M[1:], 'b', label='From forces')
M = (b.p.young*b.p.inertia)*(b.g.rho - b.g_p.rho)
plt.plot(b.g.x1_grid[1:], M[1:], 'r', label='From CST')
plt.legend()
plt.show()
plt.show()
