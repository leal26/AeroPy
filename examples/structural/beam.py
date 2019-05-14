from aeropy.geometry.parametric import poly
from aeropy.structural.stable_solution import (structure, mesh_1D, properties,
                                               boundary_conditions)
from aeropy.xfoil_module import output_reader

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pickle

abaqus_primary = pickle.load(open("save.p", "rb"), encoding='latin1')
abaqus_secondary = output_reader('secondary_variables.txt')
# sort data
abaqus_data = np.array(sorted(zip(abaqus_primary['C_U']['x'],
                                  abaqus_primary['C_U']['y'],
                                  abaqus_primary['U'][:, 0],
                                  abaqus_primary['U'][:, 1],)))
abq_x, abq_y, abq_u1, abq_u2 = abaqus_data.T
abq_y = -abq_y + .005
abq_u2 = -abq_u2
# Convert log strains into engineering strain
abaqus_secondary['LE11'] = np.exp(np.array(abaqus_secondary['LE11'])) - 1
abaqus_secondary['LE12'] = np.exp(np.array(abaqus_secondary['LE12'])) - 1
abaqus_secondary['LE22'] = np.exp(np.array(abaqus_secondary['LE22'])) - 1
coefficients = np.array([0, 0, 0, 0])

bp = properties()
bc = boundary_conditions(load=np.array([[0, -1]]))
analytical_solution = bc.concentrated_load[0][1]/(6*bp.young*bp.inertia) * \
    np.array([-1, 3, 0, 0])
mesh = mesh_1D(mesh_n=10)
curve_parent = poly(a=[0, 0, 0, 0])
curve_child = poly(a=analytical_solution)

beam = structure(curve_parent, curve_child, mesh, bp, bc)
beam.calculate_position()
strain = beam.strain()
stress = beam.stress(loading_condition='plane_stress')

# Plot beam results
plt.figure()
u = beam.u()
u1 = beam.u(diff='x1')
u2 = beam.u(diff='x2')
plt.plot(beam.r_p[0], beam.r_p[1], label='parent')
plt.scatter(beam.r_p[0], beam.r_p[1], label='parent')
plt.plot(beam.r_c[0], beam.r_c[1], label='child')
plt.scatter(beam.r_c[0], beam.r_c[1], label='child')
plt.plot(abq_x, abq_y, label='Abaqus')
plt.title('Position')
plt.grid()
plt.legend()

# Plot beam results
plt.figure()
r1_p, r1_c = beam.calculate_position(diff='x1')
r2_p, r2_c = beam.calculate_position(diff='x2')
# plt.plot(beam.r_p[0], r1_p[0], label='$r_{1,1}^p$')
plt.plot(beam.r_p[0], r1_p[1], label='$r_{2,1}^p$')
# plt.plot(beam.r_p[0], r2_p[0], label='$r_{1,2}^p$')
plt.plot(beam.r_p[0], r2_p[1], label='$r_{2,2}^p$')
# plt.plot(beam.r_p[0], r1_c[0], label='$r_{1,1}^c$')
plt.plot(beam.r_p[0], r1_c[1], label='$r_{2,1}^c$')
# plt.plot(beam.r_p[0], r2_c[0], label='$r_{1,2}^c$')
plt.plot(beam.r_p[0], r2_c[1], label='$r_{2,2}^c$')
plt.title('Position gradients')
plt.grid()
plt.legend()

# Plot beam results
plt.figure()
u = beam.u()
u1 = beam.u(diff='x1')
u2 = beam.u(diff='x2')
plt.scatter(beam.mesh.x_p, u[0], label=r'$u_1$')
plt.scatter(beam.mesh.x_p, u[1], label=r'$u_2$')
plt.plot(beam.mesh.x_p, u[0], label=r'$u_1$')
plt.plot(beam.mesh.x_p, u[1], label=r'$u_2$')
# plt.plot(abq_x, abq_u1, label=r'Abaqus $u_1$')
# plt.plot(abq_x, abq_u2, label=r'Abaqus $u_2$')
plt.title('Displacement diff')
plt.legend()


plt.figure()
plt.plot(beam.mesh.x_p, strain[0][0], label=r'$\epsilon_{11}$')
plt.plot(beam.mesh.x_p, strain[0][1], label=r'$\epsilon_{12}$')
plt.plot(beam.mesh.x_p, strain[1][1], label=r'$\epsilon_{22}$')
plt.plot(abaqus_secondary['X'], abaqus_secondary['LE11'],
         label=r'Abaqus $\epsilon_{11}$')
plt.plot(abaqus_secondary['X'], abaqus_secondary['LE12'],
         label=r'Abaqus $\epsilon_{12}$')
plt.plot(abaqus_secondary['X'], abaqus_secondary['LE22'],
         label=r'Abaqus $\epsilon_{22}$')
plt.title('Strain')
plt.legend()

plt.figure()
plt.plot(beam.mesh.x_p, stress[0][0], label=r'$\sigma_{11}$')
plt.plot(beam.mesh.x_p, stress[0][1], label=r'$\sigma_{12}$')
plt.plot(beam.mesh.x_p, stress[1][1], label=r'$\sigma_{22}$')
plt.plot(abaqus_secondary['X'], abaqus_secondary['S11'],
         label=r'Abaqus $\sigma_{11}$')
plt.plot(abaqus_secondary['X'], abaqus_secondary['S12'],
         label=r'Abaqus $\sigma_{12}$')
plt.plot(abaqus_secondary['X'], abaqus_secondary['S22'],
         label=r'Abaqus $\sigma_{22}$')
plt.legend()
plt.title('Stress')
plt.show()
