from aeropy.geometry.parametric import polynomial
from aeropy.structural.shell import shell
from aeropy.structural.stable_solution import (mesh_1D, properties,
                                               boundary_conditions,
                                               euler_bernoulle)
from optimization_tools.DOE import DOE
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# Analytical solution
bp = properties()
curve_child = polynomial([0, 0, 0, 0, 0], color = 'k')
x1 = np.linspace(0, 1, 100)
eulerBernoulle = euler_bernoulle(bp, -1, 'concentrated', curve_child, x1)
eulerBernoulle.analytical_solutions()
print(eulerBernoulle.g.D)
eulerBernoulle.g.r()
eulerBernoulle.g.plot(label='Analytical')

# Minimum Potential Energy solution
eulerBernoulle.free_energy()
eulerBernoulle.strain_energy()
eulerBernoulle.work()
eulerBernoulle.residual()
print('Analytical: ', eulerBernoulle.g.D, eulerBernoulle.R)
eulerBernoulle.minimum_potential()
print('Numerical:  ', eulerBernoulle.g.D, eulerBernoulle.R)
eulerBernoulle.g.r()
print(eulerBernoulle.g.position)
eulerBernoulle.g.plot(label='Numerical', linestyle = '--', color = '.5')
plt.legend()
plt.show()
