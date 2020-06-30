import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem

w = 1


def distributed_load(s):
    try:
        len(s)
        return w*np.ones(len(s))
    except:
        return w


def format_input(input):
    # COnsidering BC for zero derivative at the root
    return list(input) + [-input[0]]


g = CoordinateSystem.CST(D=[0, 0, 0, 0, 0], chord=1, color='b', N1=1, N2=1)
p = properties()
plt.figure()
l = loads(distributed_load=distributed_load)
s = np.linspace(0, 1, 10)

b = beam_chen(g, p, l, s)
b.parameterized_solver(format_input=format_input, x0=g.D[:-1], ignore_ends=True)

# Results from Abaqus
# abaqus_data = pickle.load(open('neutral_line.p', 'rb'))

# Euler-Bernoulle
c2 = 6*(p.length**2)*w/(24*p.young*p.inertia)
c3 = -4*(p.length)*w/(24*p.young*p.inertia)
c4 = (1)*w/(24*p.young*p.inertia)
EB_solution = - np.array([0, 0, c2, c3, c4])

curve_EB = CoordinateSystem.polynomial(D=EB_solution, chord=1, color='0.5')

curve_EB.calculate_x1(b.s)
eulerBernoulle = curve_EB.r(b.s)

[x, y] = eulerBernoulle.T

g0 = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord=1, color='k')
g0.calculate_x1(s)
g0.plot(label='Parent', zorder=0)
plt.plot(b.x, b.y, '.5', label='Child: %.3f N/m' % w, linestyle='-', lw=3, zorder=1)
plt.scatter(x, y, c='.5', label='Euler-Bernoulle: %.3f N/m' % w, edgecolors='k', zorder=2)
plt.xlim([0, 1])
# plt.scatter(abaqus_data['coord'][0:401:40,0], abaqus_data['coord'][0:401:40,1], c='g', label='FEA', edgecolors='k', zorder = 10)
plt.legend()
plt.show()
