from aeropy.geometry.parametric import poly
from aeropy.structural.stable_solution import (structure, mesh_1D, properties,
                                               boundary_conditions)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

coefficients = np.array([0, 0, 0, 0])

beam_properties = properties()
bc = boundary_conditions()
analytical_solution = bc.concentrated_load[0][0]/(beam_properties.young *
                                                  beam_properties.area)

mesh = mesh_1D(alpha=[1+analytical_solution, 1+analytical_solution],
               alpha_nodes=[0, 1], mesh_n=10)
curve_parent = poly(coefficients)
curve_child = poly(coefficients)

beam = structure(curve_parent, curve_child, mesh, beam_properties, bc)
strain = beam.strain()
stress = beam.stress()

n = 21
strain_list = analytical_solution*np.linspace(0, 2, n)
strain_matrix = np.array(list(zip(strain_list, strain_list)))
strain_x = np.linspace(0, 1, n)
strain_x = np.array(list(zip(strain_list, strain_x)))
energy, residual = beam.sweep_strains(strain_matrix, strain_x)

print(beam.find_stable())
fig = plt.figure()

# plt.yscale('log')
plt.plot(strain_list, energy, lw=2, label='Strain Energy ($U$)')
plt.plot(strain_list, residual, lw=2, label='Potential Energy ($\Pi$)')
plt.scatter(beam.opt_x, beam.opt_f)
print(beam.opt_x)
print(beam.opt_f)
plt.xlabel(r'$\epsilon_{11}$')
plt.ylabel('Energy (J)')
plt.axvline(x=analytical_solution, label='Analytical solution', c='k')
plt.legend()
plt.show()
