from aeropy.geometry.parametric import poly
from aeropy.structural.stable_solution import (structure, mesh_1D, properties,
                                               boundary_conditions)
from aeropy.xfoil_module import output_reader

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from optimization_tools.DOE import DOE
import pickle

def input_function(x):
    return np.array([0,0] + list(x))

np.set_printoptions(precision=5)

# Results from Abaqus
abaqus_data = pickle.load(open('neutral_line.p', 'rb'))

# Beam properties
bp = properties()
bc = boundary_conditions(load=np.array([[0, -1], ]))
EB_solution = bc.concentrated_load[0][1]/(6*bp.young*bp.inertia) * \
    np.array([0, 0, 3, -1])

mesh = mesh_1D(mesh_n=11, alpha=[1, 1],
               alpha_nodes=[0, 1])
curve_parent = poly(a=[0, 0, 0, 0])
curve_child = poly(a=EB_solution)

beam = structure(curve_parent, curve_child, mesh, bp, bc)
beam.calculate_position()
eulerBernoulle = beam.r_c

# DOE setup
n = 10
problem = DOE(levels=n, driver='Full Factorial')
problem.add_variable('a2', lower=-0.010, upper=-0.007, type=float)
problem.add_variable('a3', lower=0.002, upper=.005, type=float)
problem.define_points()

print(problem.domain)
coefficient_matrix = np.array(list(zip(problem.domain['a2'],
                                       problem.domain['a3'])))
# Doe itself
energy, residual = beam.sweep_geometries(coefficient_matrix,
                                         input_function = input_function,
                                         reorder = (n, n))
# Find my optimal
bounds = np.array(((-0.010,0.005),
                  (-0.005,0.005)))
print('Calculating minimum')
sol, fun = beam.find_stable(EB_solution[2:4], bounds=bounds,
                         input_type = 'Geometry',
                         loading_condition='plane_stress',
                         input_function = input_function)

fig = plt.figure()
x = np.resize(coefficient_matrix[:, 0], (n, n))
y = np.resize(coefficient_matrix[:, 1], (n, n))
print('x', x)
print('y', y)
print('residual', residual)
# residual_log = x_log_modulus_transform= np.sign(residual)*(np.log10(np.abs(residual)+1))
plt.contourf(x, y, residual, levels=100)
cbar = plt.colorbar()
cbar.set_label('Residual')
plt.scatter(sol[0], sol[1], marker = 'o', color='white', label='Energy minimum')
plt.scatter(EB_solution[2], EB_solution[3], marker = '^', color='white',
            label='Euler-bernoulli')
plt.xlabel(r'$A_2$')
plt.ylabel(r'$A_3$')

plt.legend()
plt.show()
