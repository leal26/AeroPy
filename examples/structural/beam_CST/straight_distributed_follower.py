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

x_B1 = [0, 0.111088563, 0.179019141, 0.232567106, 0.27826509, 0.317433824, 0.367033123, 0.425769313,
        0.481887716, 0.54582566, 0.599329657, 0.658042172, 0.724577611, 0.7793719, 0.86807229]
x_B0 = [0, 0.131510417, 0.2109375, 0.294270833, 0.384114583, 0.460286458,
        0.529947917, 0.598958333, 0.67578125, 0.746744792, 0.815104167, 0.889322917]
# y_B1 = [0, -0.008624495, -0.025178451, -0.043081887, -0.063603111, -0.081562353, -0.109830402, -0.143237635, -0.177948689, -0.221685023, -0.256406224, -0.298869177, -0.349063735, -0.390248233, -0.462349008]
y_B0 = [0, -0.014878806, -0.03765024, -0.068739984, -0.108779046, -0.149511719, -
        0.193474559, -0.235286058, -0.281631611, -0.328230469, -0.375198317, -0.422476963]
y_B1 = [0, -0.009918169, -0.027765799, -0.045669235, -0.064249948, -0.082856027, -0.111124076, -
        0.144531309, -0.179889201, -0.223625535, -0.257699898, -0.300809688, -0.351004246, -0.399188745, -0.465998382]
lamda = 4
# w = 400


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
w = p.young*p.inertia*lamda/p.length**3
s = np.linspace(0, 1, 21)

plt.figure()
g0 = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord=1, color='k')
g0.calculate_x1(s)
g0.plot(label='Parent', zorder=0)

l = loads(distributed_load=distributed_load, follower=False)
b = beam_chen(g, p, l, s, ignore_ends=True)
b.parameterized_solver(format_input=format_input, x0=g.D[:-1])
print('x', b.x)
print('y', b.y)
print('chord', b.g.chord)
print('deltaz', b.g.deltaz)
print('RMSE', rmse(b.x, b.y, x_B0, y_B0))
plt.plot(b.x, b.y, '.3', label=r'Child: %.3f N/m, $\beta=0$' % w, linestyle='-', lw=3, zorder=1)
plt.scatter(x_B0, y_B0, c='.3', label=r'Rao: %.3f N/m, $\beta=0$' %
            w, edgecolors='k', zorder=2, marker='s')
BREAK
print('CASE 2')
l = loads(distributed_load=distributed_load, follower=True)
b = beam_chen(g, p, l, s, ignore_ends=True)
b.parameterized_solver(format_input=format_input, x0=g.D[:-1])
plt.plot(b.x, b.y, '.5', label=r'Child: %.3f N/m, $\beta=1$' % w, linestyle='-', lw=3, zorder=1)
print('RMSE', rmse(b.x, b.y, x_B1, y_B1))
plt.scatter(x_B1, y_B1, c='.5', label=r'Rao: %.3f N/m, $\beta=1$' %
            w, edgecolors='k', zorder=2, marker='s')
plt.xlim([0, 1])
# plt.scatter(abaqus_data['coord'][0:401:40,0], abaqus_data['coord'][0:401:40,1], c='g', label='FEA', edgecolors='k', zorder = 10)
plt.legend()
plt.show()
