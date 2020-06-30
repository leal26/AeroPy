from scipy.optimize import approx_fprime
import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem

from aeropy.CST_2D import dxi_u, dxi_l, ddxi_u, ddxi_l
from aeropy.geometry.airfoil import CST


def format_input(input):
    return input[2:4]


g_fit = CoordinateSystem.CST(D=[0.00571427, 0.00428572,
                                0.00285719, -0.00571429],
                             chord=1, color='b')
g_sol = CoordinateSystem.CST(D=[5.71484605e-03, 4.28497408e-03,
                                2.85818149e-03, - 4.90351205e-08],
                             chord=1, color='b')
x = np.linspace(0.0, 1, 21)
# s = np.zeros(len(x))
# for i in range(len(x)):
#     s[i] = g.arclength(x[i])[0]

g_fit.x1_grid = x
d_fit = dxi_u(x, g_fit.D[:-1], g_fit.D[-1])
dd_fit = ddxi_u(x, g_fit.D[:-1])
d_sol = dxi_u(x, g_fit.D[:-1], g_fit.D[-1])
dd_sol = ddxi_u(x, g_sol.D[:-1])

plt.figure()
plt.plot(x, d_fit, label='dFit')
plt.plot(x, dd_fit, 'r--', label='ddFit')
plt.scatter(x, d_sol, label='dSol')
plt.scatter(x, dd_sol, label='ddSol')
plt.legend()
plt.show()
