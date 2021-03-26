import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem

from scipy.interpolate import interp1d

def rmse(x1, y1, x2, y2):
    kind = "cubic" 
    if max(x1) > max(x2):
        f = interp1d(x1, y1, kind=kind)
        predictions = y2
        targets = f(x2)
    else:
        f = interp1d(x2, y2, kind=kind)
        predictions = y1
        targets = f(x1)
    return np.sqrt(np.mean((predictions-targets)**2))
    
def format_input(input, g=None, g_p=None):
    return [0, 0] + list(input)


g = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord=1, color='b')
p = properties()

l = loads(concentrated_load=[[0, -1]], load_s=[1])
s = np.linspace(0, 1, 10)

b = beam_chen(g, p, l, s)
b.parameterized_solver(format_input, x0=b.g.D[2:])

plt.figure()
plt.plot(b.x, b.M)
plt.show()
# Results from Abaqus
abaqus_data = pickle.load(open('neutral_line.p', 'rb'))

# Euler-Bernoulle
EB_solution = l.concentrated_load[0][-1]/(6*p.young*p.inertia) * \
    np.array([0, 0, 3, -1])
print('EB', EB_solution)
print(b.g.D, b.g.zetaT)

curve_EB = CoordinateSystem.polynomial(D=EB_solution, chord=1, color='0.5')
curve_EB.calculate_x1(b.s)
eulerBernoulle = curve_EB.r(b.s)

import matplotlib
matplotlib.rcParams.update({'font.size': 14})


g0 = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord=1, color='k')
g0.calculate_x1(b.s)
g0.plot(label='Parent')
[x, y] = eulerBernoulle.T
print('RMSE', rmse(b.x, b.y, x, y))
plt.scatter(x, y, c='.5', label='Euler-Bernoulle: %.3f N' % -
            l.concentrated_load[0][1], linestyle='-', zorder=3, edgecolors='k', clip_on=False)
plt.scatter(abaqus_data['coord'][0:401:40, 0], abaqus_data['coord'][0:401:40, 1], c='.5',
            label='FEA: %.3f N' % -l.concentrated_load[0][1], edgecolors='k', zorder=4, marker="^", clip_on=False)
plt.plot(b.x, b.y, '.5', label='Child: %.3f N' % -l.concentrated_load[0][1], lw=3, zorder=2, clip_on=False)
plt.xlim([0,1])
plt.legend()


rhs = b.p.young*b.p.inertia*(b.g.rho - b.g_p.rho)
lhs = b.M
plt.figure()
print('r', min(b.r), max(b.r))
plt.plot(s, b.r, '.5', lw=3,zorder=1)
# plt.scatter(s, lhs, c='.5', edgecolors='k', zorder=20, marker="D", clip_on=False)
#plt.ylim([-b.p.young*b.p.inertia, b.p.young*b.p.inertia])
plt.xlim([0, max(s)])
plt.ylabel("Units (N.m)")
plt.xlabel("s (m)")
plt.show()

plt.show()
