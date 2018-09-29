from aeropy.geometry.parametric import poly, frame
from aeropy.structural.stable_solution import structure
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

input_type = 'z1'
x1 = np.linspace(0, 1)
x2 = 0
curve_parent = poly([0, 0, 0, 0])
curve_child = poly()
r_p = curve_parent.r(x1, x2, input_type=input_type)
r_c = curve_child.r(x1, x2, input_type=input_type)
beam = structure(curve_parent, curve_child)
u = beam.displacement(x1, input_type=input_type)
print(u)
plt.figure()
plt.plot(r_p[0], r_p[1], label='parent')
plt.plot(r_c[0], r_c[1], label='child')
Q = plt.quiver(r_p[0], r_p[1], u[0], u[1], angles='xy', scale_units='xy',
               scale=1, width=0.005, color='0.5')
plt.legend()
plt.show()

ux1 = beam.displacement(x1, input_type=input_type, diff='x1')
ux11 = beam.displacement(x1, input_type=input_type, diff='x11')
plt.figure()
plt.plot(x1, u[1], label='$u$')
plt.plot(x1, ux1[1], label=r'$\frac{\partial u}{\partial x^1}$')
plt.plot(x1, ux11[1], label=r'$\frac{\partial^2 u}{\partial (x^1)^2}$')
plt.legend()
plt.show()
