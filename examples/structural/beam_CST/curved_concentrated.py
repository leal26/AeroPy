import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


# def format_input(input, g=None, g_p=None):
#     # COnsidering BC for zero derivative at the root
#     return list(input) + [-input[0]]
def format_input(input, g=None, g_p=None):
    g.D[:-1] = input
    error = 9999
    n = g.n - 2
    dd_p = (2*n*g_p.D[-3] - 2*(g_p.N1+n)*g_p.D[-2])
    d_p = -g_p.D[-2] + g_p.D[-1]
    rho_p = dd_p/(1+d_p**2)**(3/2)
    target_length = g_p.arclength(g_p.chord)[0]
    print('target', target_length)
    while error > 1e-9:
        before = g.D[-1]
        C = rho_p*g.chord/g_p.chord
        dd = (2*n*g.D[-3] - 2*(g.N1+n)*g.D[-2])
        d = -g.D[-2] + g.D[-1]
        g.D[-1] = np.sqrt((dd/C)**(2/3)-1) + g.D[-2]
        print('terms', np.sqrt((dd/C)**(2/3)-1), g.D[-2])
        print('C', C, rho_p, g.chord, g_p.chord)
        print('dd', d, d_p, dd, dd_p)
        print('before', before)
        g.internal_variables(target_length)
        after = g.D[-1]
        print('after', after)
        error = abs(after - before)
        print('deltaxi', g.D[-1], 'error', error)
        # BREAK
    return g.D


abaqus_x = [0, 0.10398348, 0.20769666, 0.31089693, 0.4133383, 0.51480412,
            0.61510199, 0.71406609, 0.81157351, 0.90749156, 1.0017918]
abaqus_y = [0, 0.0025992838, 0.010444776, 0.023406046, 0.041411132,
            0.064284593, 0.091835514, 0.12384841, 0.16006343, 0.20029253, 0.24419737]

g = CoordinateSystem.CST(D=[-0.25, -0.25, -0.25, -0.25, 0.25], chord=1,
                         color='b', N1=1, N2=1, deltaz=0.25)
g_p = CoordinateSystem.CST(D=[-0.25, -0.25, -0.25, -0.25, 0.25], chord=1,
                           color='k', N1=1, N2=1, deltaz=0.25)

x = np.linspace(0, 1, 20)
s = np.zeros(len(x))
for i in range(len(x)):
    s[i] = g.arclength(x[i])[0]

p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[s[-1]])
# s = np.linspace(0, 1, 10)
b = beam_chen(g, p, l, s, ignore_ends=True)
b.parameterized_solver(format_input=format_input, x0=b.g.D[:-1])

g_p.calculate_x1(s)
g_p.plot(label='Parent')
plt.plot(b.x, b.y, '.5', label='Child: %.3f N' % -l.concentrated_load[0][-1], lw=3)
plt.scatter(abaqus_x, abaqus_y, c='.5', label='FEA: %.3f N' % -
            l.concentrated_load[0][-1], edgecolors='k', zorder=10, marker="^")

plt.legend()
plt.show()
