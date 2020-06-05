from aeropy.geometry.parametric import CoordinateSystem
from aeropy.structural.beam import euler_bernoulle_curvilinear
from aeropy.structural.stable_solution import (mesh_1D, properties,
                                               boundary_conditions,
                                               euler_bernoulle)
from optimization_tools.DOE import DOE
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

abaqus_x = [0, 0.050000664, 0.100002654, 0.150012404, 0.200025707, 0.250053465,
            0.300086498, 0.350139558, 0.400199026, 0.450283021, 0.50037396,
            0.550492764, 0.600618601, 0.650774479, 0.700936854, 0.751130462,
            0.801329434, 0.851559758, 0.901793659, 0.952057898, 1.002323747]
abaqus_y = [0, 0.001222659, 0.002391184, 0.006008377, 0.009577097, 0.015600263,
            0.021580685, 0.030021306, 0.038424935, 0.049294505, 0.060132839,
            0.073442832, 0.086727314, 0.102489136, 0.118231162, 0.136456192,
            0.154667079, 0.175366625, 0.196057677, 0.219242975, 0.242425412]
plt.scatter(abaqus_x, abaqus_y, c='g', label='FEA', edgecolors='k', zorder = 10)

# Analytical solution
bp = properties()
curve_child = CoordinateSystem.polynomial([0, 0, 0.25, 0.], color = 'k')
curve_parent = CoordinateSystem.polynomial([0, 0, 0.25, 0.], color ='0.5')
theta1 = np.linspace(0, curve_parent.arclength()[0], 20)
eulerBernoulle = euler_bernoulle_curvilinear(curve_parent, curve_child, bp,
                                             -1, 'concentrated', theta1)
eulerBernoulle.analytical_solutions()
eulerBernoulle.update_chord()
eulerBernoulle.g_c.calculate_x1(eulerBernoulle.theta1)
eulerBernoulle.work()
eulerBernoulle.bending_strain()
eulerBernoulle.free_energy()
eulerBernoulle.strain_energy()
eulerBernoulle.residual()
print(eulerBernoulle.g_c.D)
print(eulerBernoulle.R)
eulerBernoulle.g_c.r()
eulerBernoulle.g_p.plot(label='Baseline')
eulerBernoulle.g_c.plot(label='Analytical')

# Minimum Potential Energy solution
eulerBernoulle.free_energy()
eulerBernoulle.strain_energy()
eulerBernoulle.work()
eulerBernoulle.residual()
print('Analytical: ', eulerBernoulle.g_c.D, eulerBernoulle.R)
eulerBernoulle.minimum_potential()
print('Numerical:  ', eulerBernoulle.g_c.D, eulerBernoulle.R)
eulerBernoulle.g_c.r()
print(eulerBernoulle.g_c.position)
print('Energy', eulerBernoulle.U)
eulerBernoulle.g_c.plot(label='Numerical', linestyle = '--', color = 'b')
plt.legend()
plt.show()
