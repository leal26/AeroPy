from aeropy.geometry.parametric import polynomial
from aeropy.structural.shell import shell
from aeropy.structural.stable_solution import (mesh_1D, properties,
                                               boundary_conditions,
                                               euler_bernouille)
from optimization_tools.DOE import DOE
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

bp = properties()

chord_parent = 1

curve_child = polynomial([0, 0, 0, 0], color = 'g')

x1 = np.linspace(0, 1, 5)
eulerBernoulle = euler_bernouille(bp, 1, 'distributed', curve_child, x1)

BRAKE
# Kinematics
beam.calculate_chord()
theta_1 = np.linspace(0, 1, 10)
beam.g_p.calculate_x1(theta_1)
beam.g_c.calculate_x1(theta_1)
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
beam.g_p.curvature_tensor()
beam.g_c.curvature_tensor()
beam.g_p.r()
beam.g_c.r()
beam.calculate_change_curvature()
beam.CauchyGreen()
beam.free_energy()
beam.strain_energy()

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
#print(beam.g_p.B)
print('All energy')
print(beam.phi)
print('Bending')
print(beam.phi_B)
print('Stretching')
print(beam.phi_M)
print('change of curvature')
print(beam.rho)
print('strain energy')
print(beam.U)
plt.figure()
beam.g_p.plot(basis=True)
beam.g_c.plot(basis=True)
# plt.gca().set_aspect('equal', adjustable='box')
plt.show()

def DOE_function(inputs):
    beam.g_c.D[2] = inputs['D2']
    beam.g_c.D[3] = inputs['D3']
    beam.calculate_chord()
    beam.g_c.calculate_x1(theta_1)
    beam.g_c.basis()
    beam.g_c.metric_tensor()
    beam.calculate_strains()

    # Calculate energy
    beam.g_c.basis(diff = 'theta')
    beam.g_c.metric_tensor(diff = 'theta')
    beam.g_c.curvature_tensor()
    beam.calculate_change_curvature()
    beam.CauchyGreen()
    beam.free_energy()
    beam.strain_energy()
    beam.g_c.r()
    beam.work()
    beam.residual()
    return({'U': beam.U, 'W':beam.W, 'R':beam.R})


# Define points
print(eulerBernoulle)
problem = DOE(levels=3, driver='Full Factorial')
problem.add_variable('D2', lower=0, upper=2*eulerBernoulle[2], type=float)
problem.add_variable('D3', lower=2*eulerBernoulle[3], upper=0, type=float)
problem.define_points()

# Run for a function with dictionary as inputs
problem.run(DOE_function, tracking=True)

# Plot factor effects
# problem.plot(xlabel=['D2', 'D3'],
#              ylabel=['Strain Energy'])
# Plot domain
problem.plot_contour('D2', 'D3', 'W',  labels=['$D_2$', '$D_3$', 'Work'])
plt.scatter(eulerBernoulle[2], eulerBernoulle[3], marker = '^', color='white',
            label='Euler-bernoulli')

problem.plot_contour('D2', 'D3', 'U',  labels=['$D_2$', '$D_3$', 'Strain Energy'])
plt.scatter(eulerBernoulle[2], eulerBernoulle[3], marker = '^', color='white',
            label='Euler-bernoulli')

problem.plot_contour('D2', 'D3', 'R',  labels=['$D_2$', '$D_3$', 'Residual'])
plt.scatter(eulerBernoulle[2], eulerBernoulle[3], marker = '^', color='white',
            label='Euler-bernoulli')
plt.show()
