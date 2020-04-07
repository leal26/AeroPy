from aeropy.geometry.parametric import CoordinateSystem
from aeropy.structural.shell import shell
from aeropy.structural.stable_solution import (mesh_1D, properties,
                                               boundary_conditions,
                                               euler_bernoulle)
from optimization_tools.DOE import DOE
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pickle
import math


# bc = boundary_conditions(distributed_load=1)

# Define parent and child geometries
chord_parent = 1
load = -1.0
# Results from Abaqus

abaqus_x = [0, 0.050000664, 0.100002654, 0.150012404, 0.200025707, 0.250053465,
            0.300086498, 0.350139558, 0.400199026, 0.450283021, 0.50037396,
            0.550492764, 0.600618601, 0.650774479, 0.700936854, 0.751130462,
            0.801329434, 0.851559758, 0.901793659, 0.952057898, 1.002323747]
abaqus_y = [0, 0.001222659, 0.002391184, 0.006008377, 0.009577097, 0.015600263,
            0.021580685, 0.030021306, 0.038424935, 0.049294505, 0.060132839,
            0.073442832, 0.086727314, 0.102489136, 0.118231162, 0.136456192,
            0.154667079, 0.175366625, 0.196057677, 0.219242975, 0.242425412]

bp = properties()

curve_parent = CoordinateSystem.polynomial([0,0,.25,0], chord_parent, 'b')

bc = boundary_conditions(concentrated_load=np.array([[0, 0, load], ]),
                         load_x = [chord_parent])
curve_child = CoordinateSystem.polynomial([0,0,.25,0], chord_parent, color = '0.5')

# Define shell
beam = shell(curve_parent, curve_child, bp, bc)

# Kinematics
chord_bounds = np.array([[0, 2*chord_parent],])
beam.calculate_chord(bounds = chord_bounds)
beam.theta1 = np.linspace(0, beam.g_p.arclength()[0], 10)
beam.g_p.bounds = chord_bounds
beam.g_c.bounds = chord_bounds

beam.update_parent()
# [-7E-05, .242,0.0014]
# [0 ,.25,0]
# coefficients, results = beam.stepped_loading(x0=[.25,0],
#                                              input_function = lambda x: [0,0] + list(x),
#                                              bounds = np.array([[0.1,0.4], [0,0.02]]))
beam.minimum_potential(x0=[.25, 0], input_function = lambda x: [0,0] + list(x),
                       bounds = np.array([[0.1,0.4], [0,0.02]]))
# plt.figure()
# plt.scatter(-1*np.array(results[:,0]), -1*np.array(results[:,1]))
# plt.plot([0, 1], [0, 0.00757])
# plt.show()
print('Minimum D', beam.g_c.D, beam.R, beam.W, beam.U)
# beam.g_p.plot(label='Parent', linestyle = '-', color = 'k')
beam.g_p.plot(label='Parent')
beam.g_c.plot(label='Child', linestyle = '--')
plt.scatter(abaqus_x, abaqus_y, c='g', label='FEA', edgecolors='k', zorder = 10)
plt.legend()
plt.show()
BRAKE
def DOE_function(inputs):
    beam.g_c.D[2] = inputs['D2']
    beam.g_c.D[3] = inputs['D3']
    beam.update_child()
    return({'U': beam.U, 'W':beam.W, 'R':beam.R})


# Define points
print(eulerBernoulle)
problem = DOE(levels=20, driver='Full Factorial')
problem.add_variable('D2', lower=0, upper=2*eulerBernoulle.g.D[2], type=float)
problem.add_variable('D3', lower=2*eulerBernoulle.g.D[3], upper=0, type=float)
problem.define_points()

# Run for a function with dictionary as inputs
problem.run(DOE_function, tracking=True)

# Plot factor effects
# problem.plot(xlabel=['D2', 'D3'],
#              ylabel=['Strain Energy'])
# Plot domain
eulerBernoulle.analytical_solutions()
print(eulerBernoulle.g.D)
problem.plot_contour('D2', 'D3', 'W',  labels=['$D_2$', '$D_3$', 'Work'])
plt.scatter(eulerBernoulle.g.D[2], eulerBernoulle.g.D[3], marker = '^',
            color='white', label='Euler-bernoulli')
plt.scatter(x[0], x[1], marker = 's', color='k', label='Minimum Potential')
plt.legend()

problem.plot_contour('D2', 'D3', 'U',  labels=['$D_2$', '$D_3$', 'Strain Energy'])
plt.scatter(eulerBernoulle.g.D[2], eulerBernoulle.g.D[3], marker = '^',
            color='white', label='Euler-bernoulli')
plt.scatter(x[0], x[1], marker = 's', color='k', label='Minimum Potential')
plt.legend()

problem.plot_contour('D2', 'D3', 'R',  labels=['$D_2$', '$D_3$', 'Residual'])
plt.scatter(eulerBernoulle.g.D[2], eulerBernoulle.g.D[3], marker = '^',
            color='white', label='Euler-bernoulli')
plt.scatter(x[0], x[1], marker = 's', color='k', label='Minimum Potential')
plt.legend()
plt.show()
