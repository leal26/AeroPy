from aeropy.geometry.parametric import polynomial
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
bc = boundary_conditions(load=np.array([[0, -1], ]))
EB_solution = bc.concentrated_load[0][1]/(6*bp.young*bp.inertia) * \
    np.array([0, 0, 3, -1])

curve_parent = polynomial(a=[0, 0, 0, 0])
curve_child = polynomial(a=EB_solution)

beam = shell(curve_parent, curve_child, bp, bc)
beam.calculate_position()
eulerBernoulle = beam.r_c

# Find stable solution
bounds = np.array(((-0.02,0.02),
                  (-0.02,0.02)))

x,fun = beam.find_stable(beam.g_c.a[2:4], bounds=bounds, input_type = 'Geometry',
                 loading_condition='plane_stress',
                 input_function = input_function)
beam.g_c.a = input_function(x)
beam.mesh.mesh_child()
beam.calculate_position()
beam.strain()
beam.stress(loading_condition='plane_stress')

print('x', x)
print('f', beam.opt_f)
print('Residual', beam.residual(x, input_type = 'Geometry',
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
