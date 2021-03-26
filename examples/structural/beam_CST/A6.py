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
    # COnsidering BC for zero derivative at the root
    return list(input) + [-input[0]]

w = 400

def distributed_load(s):
    try:
        len(s)
        return w*np.ones(len(s))
    except:
        return w


abaqus = np.transpose(np.loadtxt("A6_neutral_line.csv", delimiter=","))
abaqus = abaqus[abaqus[:,1].argsort()]
  
# Ghuku paper
distributed_load = distributed_load
chord_parent = 0.4385
width = 0.0385
thickness = 0.00625
clamped_load = 20.9634
experiment = {0: [0., 0., 0.5492, 0.2000],
              138.321+clamped_load: [0., 0., 0.3746, 0.1753],
              211.896+clamped_load: [0., 0., 0.2475, 0.2864], }
# 287.433: [0., 0., 0.1287, 0.3582],
# 460.089: [0., 0., -0.1891, 0.6087],}
p = properties(young=200e9, dimensions=[width, thickness])

exp_F = [0.000, 20.9634, 83.7474, 159.2844, 232.8594, 308.3964, 382.9542, 459.4704, 481.0524]
exp_delta = [0, 0.0047, 0.0135, 0.03098, 0.0451, 0.0622, 0.0793, 0.0977, 0.103]
g0 = CoordinateSystem.polynomial(D=experiment[list(experiment.keys())[
    0]], chord=chord_parent, color='b')

colors = ['0.0', '0.3', '0.5', '0.7']
load_keys = list(experiment.keys())

x = np.linspace(0, chord_parent, 10)
s = np.zeros(len(x))
for i in range(len(x)):
    s[i] = g0.arclength(x[i])[0]
s[0] = 0

g = CoordinateSystem.CST(D=[-0.27928065, - 0.29209946, - 0.30491828, - 0.3177371,
                            0.27928065], chord=chord_parent, color='b', deltaz=0.27928065*chord_parent, N1=1.0, N2=1.0)

for i in [0,1]:
    print('Load: ', load_keys[i])
    load = load_keys[i]
    l = loads(concentrated_load=[[0, -load]], load_s=[s[-1]], distributed_load=distributed_load)
    b = beam_chen(g, p, l, s, ignore_ends=True)
    # b.iterative_solver()
    b.parameterized_solver(format_input=format_input, x0=b.g.D[:-1])

    if i == 0:
        plt.plot(b.x, b.y, colors[i], label='Parent', linestyle='-',
                 lw=3, zorder=1)
        gg = CoordinateSystem.polynomial(D=experiment[load_keys[i]], chord=chord_parent)
        gg.calculate_x1(s)
        gg.plot(label='Experiment: %.3f N' %
                load_keys[i], color=colors[i], scatter=True, marker="D")
    if i > 0:
        plt.plot(b.x, b.y, colors[i], label='Child: %.3f N' % load, linestyle='-',
                 lw=3, zorder=1)
        plt.scatter(abaqus[:, 0][0::4], abaqus[:, 1][0::4], c='.5', edgecolors='k', zorder=10, marker="^")
        # print('RMSE', rmse(b.x, b.y, abaqus[:, 0], abaqus[:, 1]))
plt.legend()
plt.xlim([0, max(b.x)])
plt.show()
