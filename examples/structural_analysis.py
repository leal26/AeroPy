from aeropy.geometry.parametric import poly, frame
from aeropy.structural.stable_solution import structure
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

input_type = 'x1'
x1 = np.linspace(0, 1)
x2 = 0
b = 0.01
h = 0.01
P = 1000
L = 1
E = 70e9
Area = h*b
Inertia = b*h**3/12.

# For Euler Bernouille
# a = P/(6*E*Inertia)*np.array([-1, 3*L, 0, 0])
# For axial deformation
a = np.array([0, 0, 0, 0])
scaling = 1+1  # P/Area/E

curve_parent = poly([0, 0, 0, 0])
curve_child = poly(a, scaling)
r_p = curve_parent.r(x1, x2, input_type=input_type)
r_c = curve_child.r(x1, x2, input_type=input_type)
beam = structure(curve_parent, curve_child, area=b*h, young=E)
u = beam.u(x1, input_type=input_type)
strain = beam.strain(x1)
stress = beam.stress()
print(stress)
plt.figure()
plt.plot(r_p[0], r_p[1], label='parent', color='r')
plt.plot(r_c[0], r_c[1], label='child', color='b')
Q = plt.quiver(r_p[0], r_p[1], u[0], u[1], angles='xy', scale_units='xy',
               scale=1, width=0.005, color='0.5')
plt.legend()
plt.show()

ux1 = beam.u(r_p[0], input_type=input_type, diff='x1')
# ux11 = beam.u(r_p[0], input_type=input_type, diff='x11')

plt.figure()
plt.plot(r_p[0], u[0], label='$u^1$')
plt.plot(r_p[0], u[1], label='$u^2$')
plt.plot(r_p[0], ux1[0], label=r'$\frac{\partial u^1}{\partial x^1}$')
plt.plot(r_p[0], ux1[1], label=r'$\frac{\partial u^2}{\partial x^1}$')
# plt.plot(r_p[0], ux11[1], label=r'$\frac{\partial^2 u}{\partial (x^1)^2}$')
plt.legend()
plt.show()

plt.figure()
plt.plot(x1, stress[0][0], label=r'$\sigma_{calculated}$')

plt.plot(x1, -np.ones(len(x1))*P/Area, label=r'$\sigma_{EB}$')
# /2/Inertia*h**2/4.
plt.legend()
plt.show()
