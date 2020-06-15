from aeropy.geometry.parametric import CoordinateSystem
from aeropy.structural.shell import shell
from aeropy.structural.stable_solution import (mesh_1D, properties,
                                               boundary_conditions,
                                               euler_bernoulle)
from optimization_tools.DOE import DOE
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

bp = properties()
bc = boundary_conditions(concentrated_load=np.array([[0, 0, -1], ]))
# bc = boundary_conditions(distributed_load=1)

# Define parent and child geometries
straight = np.array([0, 0, 0, 0, 0])
chord_parent = 1

curve_parent = CoordinateSystem.polynomial(np.copy(straight), chord_parent, 'b')
curve_child = CoordinateSystem.polynomial(np.array([0, 0, 0, 0, 0]), chord_parent, color = 'g')
eulerBernoulle = euler_bernoulle(bp, bc.concentrated_load[0][2], 'concentrated',
                                 curve_child)
eulerBernoulle.analytical_solutions()

# Define shell
beam = shell(curve_parent, eulerBernoulle.g, bp, bc)

# Kinematics
beam.calculate_chord()
beam.theta1 = np.linspace(0, 1, 10)
beam.g_p.calculate_x1(beam.theta1)
beam.g_c.calculate_x1(beam.theta1)
beam.g_p.basis()
beam.g_c.basis()
beam.g_p.metric_tensor()
beam.g_c.metric_tensor()
beam.calculate_strains()
# Calculate energy
beam.g_p.basis(diff = 'theta')
beam.g_c.basis(diff = 'theta')
beam.g_p.metric_tensor(diff = 'theta')
beam.g_c.metric_tensor(diff = 'theta')
# print('parent')
# print(beam.g_p.A)
# print(beam.g_p.dA)
# print('child')
# print('a')
# print(beam.g_c.a)
# print('A')
# print(beam.g_c.A)
# print(beam.g_c.dA)
beam.g_p.curvature_tensor()
beam.g_c.curvature_tensor()
beam.g_p.r()
beam.g_c.r()


beam.calculate_change_curvature()
beam.CauchyGreen()
beam.free_energy()
beam.strain_energy()
beam.work()
beam.residual()

# print(beam.g_c.D)

# print('Parent vectors')
# print(beam.g_p.A)
# print('Child vector')
# print(beam.g_c.A)
# print('C')
# print(beam.C)
# print('rho')
# print(beam.rho)
# print('Child B')
# print(beam.g_c.B)
# print('Parent B')
# print(beam.g_p.B)
# print('All energy')
# print(beam.phi)
# print('Bending')
# print(beam.phi_B)
# print('Stretching')
# print(beam.phi_M)
# print('change of curvature')
# print(beam.rho)
# print('strain')
# print(beam.gamma)
print('strain energy')
print(beam.U)
print(beam.W)
print(beam.R)
x, R = beam.minimum_potential()
print(x, R)
plt.figure()
beam.g_p.plot(label='Parent')
beam.g_c.plot(label='Euler-Bernoulle', color = 'k')
beam.g_c.D[2:] = x
beam.update_child()
beam.g_c.plot(label='Minimum Potential', color = '.5', linestyle= '--')
plt.legend()
plt.show()
# BRAKE
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
