from aeropy.geometry.parametric import poly
from aeropy.structural.stable_solution import (structure, mesh_1D, properties,
                                               boundary_conditions)
from aeropy.xfoil_module import output_reader

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from optimization_tools.DOE import DOE
import pickle

def input_function(x):
    return np.array([0] + list(x))

np.set_printoptions(precision=4)

# Results from Abaqus
abaqus_primary = pickle.load(open('save.p', 'rb'), encoding='latin1')
abaqus_secondary = output_reader('secondary_variables.txt')

# sort data
abaqus_data = np.array(sorted(zip(abaqus_primary['C_U']['x'],
                                  abaqus_primary['C_U']['y'],
                                  abaqus_primary['U'][:, 0],
                                  abaqus_primary['U'][:, 1],)))
abq_x, abq_y, abq_u1, abq_u2 = abaqus_data.T
abq_y = -abq_y + .005
abq_u2 = -abq_u2

# Beam properties
bp = properties()
bc = boundary_conditions(load=np.array([[0, -1], ]))
analytical_solution = bc.concentrated_load[0][1]/(6*bp.young*bp.inertia) * \
    np.array([0, 0, 3, -1])

mesh = mesh_1D(mesh_n=10, alpha=[1, 1],
               alpha_nodes=[0, 1])
curve_parent = poly(a=[0, 0, 0, 0])
curve_child = poly(a=analytical_solution)

beam = structure(curve_parent, curve_child, mesh, bp, bc)
beam.calculate_position()
eulerBernoulle = beam.r_c
strain = beam.strain()
stress = beam.stress(loading_condition='plane_stress')

# Find stable solution
bounds = np.array(((-0.1,0.1),
                  (-0.1,0.1),
                  (-0.1,0.1)))

x,fun = beam.find_stable(beam.g_c.a[1:4], bounds=bounds, input_type = 'Geometry',
                 loading_condition='plane_stress',
                 input_function = input_function)
beam.calculate_position()
print('x', x)
print(beam.opt_f)
print(beam.residual(beam.g_c.a[1:4], input_type = 'Geometry',
                     loading_condition = 'plane_stress',
                     input_function = input_function))
# Plot beam results
plt.figure()
u = beam.u()
u1 = beam.u(diff='x1')
u2 = beam.u(diff='x2')
plt.plot(beam.r_p[0], beam.r_p[1], label='parent')
plt.scatter(beam.r_p[0], beam.r_p[1], label='parent')
plt.plot(beam.r_c[0], beam.r_c[1], label='child')
plt.scatter(beam.r_c[0], beam.r_c[1], label='child')
plt.plot(eulerBernoulle[0], eulerBernoulle[1], label='Euler-Bernoulle')
plt.scatter(eulerBernoulle[0], eulerBernoulle[1], label='Euler-Bernoulle')
plt.plot(abq_x, abq_y, label='Abaqus')
plt.title('Position')
plt.grid()
plt.legend()

plt.figure()
plt.plot(beam.opt_f)
plt.show()


BRAKE
n = 10

problem = DOE(levels=n, driver='Full Factorial')
problem.add_variable('a1', lower=0, upper=.001, type=float)
problem.define_points()
print(problem.domain)
strain_matrix = np.array(list(zip(problem.domain['a1'],
                                  problem.domain['a1'])))
strain_x = np.array(list(zip(0*np.ones(n),
                             1*np.ones(n))))
energy, residual = beam.sweep_geometries(strain_matrix, strain_x)

beam.find_stable()

fig = plt.figure()

# plt.yscale('log')
plt.plot(strain_list, energy, lw=2, label='Strain Energy ($U$)')
plt.plot(strain_list, residual, lw=2, label='Potential Energy ($\Pi$)')
plt.scatter(beam.opt_x, beam.opt_f)
plt.xlabel(r'$\epsilon_{11}$')
plt.ylabel('Energy (J)')
plt.axvline(x=analytical_solution, label='Analytical solution', c='k')
plt.legend()
plt.show()

# fig = plt.figure()
# print('strain', strain_matrix)
# x = np.resize(strain_matrix[:, 0], (n, n))
# y = np.resize(strain_matrix[:, 1], (n, n))
# print('x', x)
# print('y', y)
# print('residual', residual)
# plt.contourf(x, y, residual)
# plt.xlabel(r'$\epsilon_a$')
# plt.ylabel(r'$\epsilon_b$')
# cbar = plt.colorbar()
# cbar.set_label('Residual')
# plt.show()
