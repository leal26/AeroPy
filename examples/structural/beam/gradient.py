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
                             chord=1, color='b', N1=1, N2=1)
g_sol = CoordinateSystem.CST(D=[5.71484605e-03, 4.28497408e-03,
                                2.85818149e-03, - 4.90351205e-08],
                             chord=1, color='b', N1=1, N2=1)

s = np.linspace(0, 1, 10)
p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[1])
b_fit = beam_chen(g_fit, p, l, s, ignore_ends=True)
b_fit.g.internal_variables(b_fit.length)
b_fit.g.calculate_x1(b_fit.s)
b_fit.integral_ends()
b_fit.x = b_fit.g.x1_grid
b_fit.y = b_fit.g.x3(b_fit.x)

b_sol = beam_chen(g_sol, p, l, s, ignore_ends=True)
b_sol.g.internal_variables(b_sol.length)
b_sol.g.calculate_x1(b_fit.s)
b_sol.integral_ends()
b_sol.x = b_sol.g.x1_grid
b_sol.y = b_sol.g.x3(b_fit.x)

x = b_fit.g.x1_grid
print(x)
d_fit = b_fit.g.x3(x, diff='x1')
print(d_fit)
dd_fit = b_fit.g.x3(x, diff='x11')
d_sol = b_sol.g.x3(x, diff='x1')
dd_sol = b_sol.g.x3(x, diff='x11')

plt.figure()
plt.plot(x, d_fit, label='dFit')
plt.plot(x, dd_fit, 'r--', label='ddFit')
plt.scatter(x, d_sol, label='dSol')
plt.scatter(x, dd_sol, label='ddSol')
plt.legend()
plt.show()
