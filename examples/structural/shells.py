from aeropy.geometry.parametric import polynomial
from aeropy.structural.shell import shell
from aeropy.structural.stable_solution import (mesh_1D, properties,
                                               boundary_conditions)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

bp = properties()
bc = boundary_conditions(load=np.array([[0, -1], ]))
analytical_solution = bc.concentrated_load[0][0]/(bp.young*bp.area)

# Define parent and child geometries
straight = np.array([0, 0, 0, 0])
eulerBernoulle = bc.concentrated_load[0][1]/(6*bp.young*bp.inertia) * \
    np.array([0, 0, 3, -1])
chord_parent = 1

curve_parent = polynomial(straight, chord_parent, color = 'b')
curve_child = polynomial(eulerBernoulle, color = 'g')

# Define shell
beam = shell(curve_parent, curve_child, bp, bc)

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
print('A')
print(beam.g_p.A)
print('a')
print(beam.g_c.A)
print('strains')
print(beam.eps)

plt.figure()
beam.g_p.plot(basis=True)
beam.g_c.plot(basis=True)
# plt.gca().set_aspect('equal', adjustable='box')
plt.show()
