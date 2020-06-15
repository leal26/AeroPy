from aeropy.geometry.parametric import poly
from aeropy.structural.stable_solution import (structure, mesh_1D, properties,
                                               boundary_conditions)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from optimization_tools.DOE import DOE

coefficients = np.array([0, 0, 0, 0])

beam_properties = properties()
bc = boundary_conditions(load=np.array([[0, 0], ]))
analytical_solution = bc.concentrated_load[0][0]/(beam_properties.young *
                                                  beam_properties.area)

mesh = mesh_1D(alpha=[1+analytical_solution, 1+analytical_solution, 1+analytical_solution],
               alpha_nodes=[0, .5, 1], mesh_n=10)
curve_parent = poly(coefficients)
curve_child = poly(coefficients)

beam = structure(curve_parent, curve_child, mesh, beam_properties, bc)
strain = beam.strain()
stress = beam.stress()

n = 10

problem = DOE(levels=n, driver='Full Factorial')
problem.add_variable('a1', lower=0, upper=.001, type=float)
problem.add_variable('a2', lower=0, upper=.001, type=float)
problem.define_points()
print(problem.domain)
strain_matrix = np.array(list(zip(problem.domain['a1'],
                                  problem.domain['a2'],
                                  problem.domain['a2'])))
strain_x = np.array(list(zip(0*np.ones(n*n),
                             0.5*np.ones(n*n),
                             1*np.ones(n*n))))
energy, residual = beam.sweep_strains(strain_matrix, strain_x, (n, n))

fig = plt.figure()
print('strain', strain_matrix)
x = np.resize(strain_matrix[:, 0], (n, n))
y = np.resize(strain_matrix[:, 1], (n, n))
print('x', x)
print('y', y)
print('residual', residual)
plt.contourf(x, y, residual)
plt.xlabel(r'$\epsilon_a$')
plt.ylabel(r'$\epsilon_b$')
cbar = plt.colorbar()
cbar.set_label('Residual')
plt.show()
