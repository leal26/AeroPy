from aeropy.geometry.parametric import poly, frame
from aeropy.structural.stable_solution import structure
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

x1 = np.linspace(0, 1, 10)
x2 = 0
b = 0.01
h = 0.01
P = 10000
L = 1
E = 70e9
Area = h*b
Inertia = b*h**3/12.

# For Euler Bernouille
# a = P/(6*E*Inertia)*np.array([-1, 3*L, 0, 0])
# For axial deformation
a = np.array([0, 0, 0, 0])
alpha = [1+P/Area/E]

curve_parent = poly(a, config='parent')
curve_child = poly(a, alpha, config='child')
r_p = curve_parent.r(x1, x2)
r_c = curve_child.r(x1, x2)
beam = structure(curve_parent, curve_child, area=Area, young=E, load=P)
u = beam.u(x1)
strain = beam.strain(x1)
stress = beam.stress()
# print('strain', strain)
# print('target', beam.load/beam.area/beam.young)
#
# plt.figure()
# plt.plot(r_p[0], r_p[1], label='parent', color='r')
# plt.plot(r_c[0], r_c[1], label='child', color='b')
# Q = plt.quiver(r_p[0], r_p[1], u[0], u[1], angles='xy', scale_units='xy',
#                scale=1, width=0.005, color='0.5')
# plt.legend()
# plt.show()
#
# ux1 = beam.u(r_p[0], diff='x1')
#
# plt.figure()
# plt.plot(x1, u[0], label='$u^1$')
# plt.plot(x1, u[1], label='$u^2$')
# plt.plot(x1, ux1[0], label=r'$\frac{\partial u^1}{\partial x^1}$')
# plt.plot(x1, ux1[1], label=r'$\frac{\partial u^2}{\partial x^1}$')
# plt.plot(x1, np.ones(len(x1))*P/Area/E, label='target')
# # plt.plot(r_p[0], ux11[1], label=r'$\frac{\partial^2 u}{\partial (x^1)^2}$')
# plt.legend()
# plt.show()
#
# plt.figure()
# plt.plot(x1, strain[0][0], label=r'$\epsilon_{calculated}$')
# print('strain target', P/Area/E)
# plt.plot(x1, np.ones(len(x1))*P/Area/E, label='target')
# plt.ylim(0, 1.1*max(strain[0][0]))
# plt.legend()
# plt.show()
#
# plt.figure()
# plt.plot(x1, stress[0][0], label=r'$\sigma_{calculated}$')
# plt.plot(x1, np.ones(len(x1))*P/Area, label=r'$\sigma_{EB}$')
# plt.ylim(0, 1.1*max(stress[0][0]))
# plt.legend()
# plt.show()

# [stable_strain, stable_residual] = beam.find_stable()
# print(stable_strain, stable_residual)
alpha_list = 1+P/Area/E*np.linspace(0, 2, 51)
energy_list = []
residual_list = []
for alpha in alpha_list:
    beam.g_c.alpha = alpha
    beam.strain(x1)
    beam.stress()
    energy_list.append(beam.strain_energy())
    residual_list.append(beam.residual())


fig = plt.figure()

# plt.yscale('log')
plt.plot(alpha_list-1, energy_list, lw=2, label='Strain energy')
plt.plot(alpha_list-1, residual_list, lw=2, label='Residual')
# plt.scatter(stable_strain, stable_residual)
plt.xlabel(r'$\epsilon_{11}$')
plt.ylabel('Energy')
plt.axvline(x=P/Area/E, label='Expected strain', c='k')
plt.legend()
plt.show()
