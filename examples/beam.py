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

fig = plt.figure()
plt.plot(mesh.x_p, beam.u()[0])
plt.plot(mesh.x_p, beam.u()[1])
print(beam.epsilon)

strain_list = analytical_solution*np.linspace(0, 2, 21)
energy, residual = beam.sweep_geometries(strain_list)

fig = plt.figure()

print('strains', strain_list)
print('energy', energy)
print('residual', residual)
# plt.yscale('log')
plt.plot(strain_list, energy, lw=2, label='Strain energy')
plt.plot(strain_list, residual, lw=2, label='Residual')
# plt.scatter(stable_strain, stable_residual)
plt.xlabel(r'$\epsilon_{11}$')
plt.ylabel('Energy')
plt.axvline(x=analytical_solution, label='Expected strain', c='k')
plt.legend()
plt.show()
