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

bp = properties()
bc = boundary_conditions(concentrated_load=np.array([[0, 0, -1], ]))
# bc = boundary_conditions(distributed_load=1)

abaqus_data = pickle.load(open('curved_neutral_line.p', 'rb'))

# Define parent and child geometries
chord_parent = math.sqrt(3)

curve_parent = CoordinateSystem.cylindrical([2], chord = chord_parent, color = 'b')
curve_child = CoordinateSystem.cylindrical([2], chord = chord_parent, color = 'g')

# Define shell
beam = shell(curve_parent, curve_child, bp, bc)

# Setting mesh and bounds
bounds =  [[.001, beam.g_c.D[0]],]
beam.calculate_chord(bounds = bounds)
beam.theta1 = np.linspace(0, beam.g_p.arclength()[0], 10)
beam.g_p.bounds = bounds
beam.g_c.bounds = bounds

beam.update_parent()
x, R = beam.minimum_potential(x0 = beam.g_p.D, bounds = [[.5*beam.g_p.D[0], 2*beam.g_p.D[0]],])

print(x, R)
plt.figure()
beam.g_p.plot(label='Parent')
# beam.g_c.plot(label='Euler-Bernoulle', color = 'k')
# beam.g_c.D[2:] = x
beam.update_child()
beam.g_c.plot(label='Minimum Potential', color = '.5', linestyle= '-')
plt.plot(abaqus_data['coord'][:,0], abaqus_data['coord'][:,1], 'r', label='FEA', lw=3)
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
