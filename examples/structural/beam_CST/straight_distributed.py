import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


w = 160

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


def distributed_load(s):
    try:
        len(s)
        return w*np.ones(len(s))
    except:
        return w


def format_input(input, g=None, g_p=None):
    # COnsidering BC for zero derivative at the root
    return list(input) + [-input[0]]


g = CoordinateSystem.CST(D=[0, 0, 0, 0, 0], chord=1, color='b', N1=1, N2=1)
p = properties()
plt.figure()
l = loads(distributed_load=distributed_load, follower=True)
s = np.linspace(0, 1, 41)

b = beam_chen(g, p, l, s, ignore_ends=False)
b.parameterized_solver(format_input=format_input, x0=g.D[:-1])

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
abaqus = np.loadtxt('straight_distributed.csv', delimiter=',')
plt.scatter(abaqus[0, :], abaqus[1, :], c='.5', label='Euler-Bernoulle: %.3f N/m' % w, edgecolors='k', zorder=2)
# plt.scatter(x, y, c='.5', label='Euler-Bernoulle: %.3f N/m' % w, edgecolors='k', zorder=2)
plt.xlim([0, 1])
# plt.scatter(abaqus_data['coord'][0:401:40,0], abaqus_data['coord'][0:401:40,1], c='g', label='FEA', edgecolors='k', zorder = 10)
plt.legend()

print('x', b.x)
print('s', b.s)
print('s0', b.s0)
rhs = b.p.young*b.p.inertia*(b.g.rho - b.g_p.rho)
lhs = b.M
plt.figure()
plt.plot(s, rhs, '.5', lw=3,zorder=1, clip_on=False)
plt.scatter(s, lhs, c='.5', edgecolors='k', zorder=20, marker="D", clip_on=False)
#plt.ylim([-b.p.young*b.p.inertia, b.p.young*b.p.inertia])
plt.xlim([0, max(s)])
plt.ylabel("Units (N.m)")
plt.xlabel("s (m)")
plt.show()
