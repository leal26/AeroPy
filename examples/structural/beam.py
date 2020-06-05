from aeropy.geometry.parametric import CoordinateSystem
from aeropy.structural.stable_solution import (properties, boundary_conditions)
from aeropy.structural.shell import shell
from aeropy.xfoil_module import output_reader

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
bc = boundary_conditions(concentrated_load=np.array([[0.0, 0.0, -1.0], ]))
EB_solution = bc.concentrated_load[0][2]/(6*bp.young*bp.inertia) * \
    np.array([0, 0, 3, -1])

curve_parent = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord = 1, color = 'b')
curve_child = CoordinateSystem.polynomial(D=EB_solution, chord = 1, color  ='0.5')

beam = shell(curve_parent, curve_child, bp, bc)
chord_bounds = [[0., 2],]
beam.g_p.bounds = chord_bounds
beam.g_c.bounds = chord_bounds
beam.theta1 = np.linspace(0, beam.g_p.arclength()[0], 20)
beam.update_parent()
beam.update_child()
eulerBernoulle = beam.g_c.r(beam.g_c.x1_grid)

# Find stable solution
bounds = np.array(((-0.02,0.02),
                  (-0.02,0.02)))

print(EB_solution[2:4])
beam.minimum_potential(x0=[0,0], input_function = lambda x: [0,0] + list(x),
                       bounds = bounds)
# coefficients, results, arc_length = beam.stepped_loading(x0=[0,0],
#                                              input_function = lambda x: [0,0] + list(x),
#                                              bounds = bounds)
# print('arc_length')
# print(arc_length)
# plt.figure()
# plt.scatter(-1*np.array(results[:,0]), -1*np.array(results[:,1]))
# plt.plot([0, 1], [0, 0.00564])
# plt.show()
#
# plt.figure()
# plt.plot(arc_length[:,0], arc_length[:,1])
# plt.show()

# beam.g_c.a = input_function(x)
# beam.mesh.mesh_child()
# beam.calculate_position()
# beam.strain()
# beam.stress(loading_condition='plane_stress')
#
# print('x', x)
# print('f', beam.opt_f)
# print('Residual', beam.residual(x, input_type = 'Geometry',
#                                 loading_condition = 'plane_stress',
#                                 input_function = input_function))
print('Strain energy for min: ', beam.U)
[x,y,z] = eulerBernoulle.T
beam.g_p.plot(label='Parent')
plt.plot(x,z, 'k', lw = 3, label='Euler-Bernoulle', linestyle = '-', zorder=0)
beam.g_c.plot(label='Child', linestyle = '--')



plt.scatter(abaqus_data['coord'][0:401:40,0], abaqus_data['coord'][0:401:40,1], c='g', label='FEA', edgecolors='k', zorder = 10)
plt.legend()
plt.show()

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
plt.plot(abaqus_data['coord'][:,0], abaqus_data['coord'][:,1], label='Abaqus')
plt.title('Position')
plt.grid()
plt.legend()

plt.figure()
plt.plot(beam.opt_f)
plt.ylabel('Energy (J)')
plt.xlabel('Iteration')


plt.figure()
plt.plot(beam.mesh.x_p, beam.epsilon[0][0], label=r'$\epsilon_{11}$')
plt.plot(beam.mesh.x_p, beam.epsilon[0][1], label=r'$\epsilon_{12}$')
plt.plot(beam.mesh.x_p, beam.epsilon[1][1], label=r'$\epsilon_{22}$')
plt.plot(abaqus_data['coord'][:,0], abaqus_data['epsilon'][0][0],
         label=r'Abaqus $\epsilon_{11}$')
plt.plot(abaqus_data['coord'][:,0], abaqus_data['epsilon'][0][1],
         label=r'Abaqus $\epsilon_{12}$')
plt.plot(abaqus_data['coord'][:,0], abaqus_data['epsilon'][1][1],
         label=r'Abaqus $\epsilon_{22}$')
plt.title('Strain')
plt.legend()

plt.figure()
plt.plot(beam.mesh.x_p, beam.sigma[0][0], label=r'$\sigma_{11}$')
plt.plot(beam.mesh.x_p, beam.sigma[0][1], label=r'$\sigma_{12}$')
plt.plot(beam.mesh.x_p, beam.sigma[1][1], label=r'$\sigma_{22}$')
plt.plot(abaqus_data['coord'][:,0], abaqus_data['sigma'][0][0],
         label=r'Abaqus $\sigma_{11}$')
plt.plot(abaqus_data['coord'][:,0], abaqus_data['sigma'][0][0],
         label=r'Abaqus $\sigma_{12}$')
plt.plot(abaqus_data['coord'][:,0], abaqus_data['sigma'][0][0],
         label=r'Abaqus $\sigma_{22}$')
plt.legend()
plt.title('Stress')
plt.show()
