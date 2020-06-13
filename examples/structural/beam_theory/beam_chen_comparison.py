import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def format_input(input):
    return input[2:4]
# Ghuku paper
# distributed_load = None
# chord_parent = 0.4385
# width = 0.0385
# thickness = 0.00625
# clamped_load = 20.9634
# experiment = {0: [0.,0.,0.5492,0.2000],
#               138.321+clamped_load: [0., 0., 0.3746, 0.1753],
#               211.896+clamped_load: [0., 0., 0.2475, 0.2864],}
#               # 287.433: [0., 0., 0.1287, 0.3582],
#               # 460.089: [0., 0., -0.1891, 0.6087],}
# p = properties(young=200e9, dimensions=[width, thickness])
# exp_F = [0.000, 20.9634, 83.7474, 159.2844, 232.8594, 308.3964, 382.9542, 459.4704, 481.0524]
# exp_delta = [0, 0.0047, 0.0135, 0.03098, 0.0451, 0.0622, 0.0793, 0.0977, 0.103]
# g = CoordinateSystem.polynomial(D=experiment[list(experiment.keys())[0]], chord = chord_parent, color = 'b')

# Chen
def distributed_load(s):
    try:
        len(s)
        return 0.758*np.ones(len(s))
    except:
        return 0.758

chord_parent = 0.4
width = 0.025
thickness = 0.0004
D_0 = [0, 0, 0, 0]

experiment = {0: [0, 0,-1.05540301,1.19321522],
              0.098: [0, 0, -1.72387251, 1.5985958],
              0.196: [0, 0, -2.3139154, 1.77815251],}

p = properties(young=194.3e9, dimensions=[width, thickness])
exp_F = [0.000, 0.098, 0.196, 0.294, 0.392, 0.490, 0.588]
exp_delta = [0.089, 0.149, 0.195, 0.227, 0.251, 0.268, 0.281]
g = CoordinateSystem.polynomial(D=D_0, chord = chord_parent, color = 'b')

colors = ['0.0','0.3', '0.5', '0.7']
load_keys = list(experiment.keys())

x = np.linspace(0, chord_parent, 10)
s = np.zeros(len(x))
for i in range(len(x)):
    s[i] = g.arclength(x[i])[0]

for i in range(len(load_keys)):
    print('Load: ', load_keys[i])
    load = load_keys[i]
    l = loads(concentrated_load = [[0, -load]], load_s = [s[-1]], distributed_load = distributed_load)
    b = beam_chen(g, p, l, s)
    # b.iterative_solver()
    b.parameterized_solver(format_input = format_input)
    plt.plot(b.x, b.y, colors[i], label='Model: %.3f N' % load, linestyle = '--',
             lw = 3)

    gg = CoordinateSystem.polynomial(D=experiment[load_keys[i]], chord = chord_parent, color = 'b')
    gg.calculate_x1(s)
    gg.plot(label='Experiment: %.3f N' % load_keys[i], color = colors[i], scatter = True)
plt.legend()
plt.show()
