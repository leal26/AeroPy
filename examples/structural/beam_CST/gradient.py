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


g_fit = CoordinateSystem.CST(D=[0.11621803608468839, 0.11164077707202186, 0.08581012060178156, 0.11474758149621056, -0.005855939703725668],
                             chord=1, color='b', N1=.5, N2=1, tol=0.0014553076272791395)
g_sol = CoordinateSystem.CST(D=[0.11314995843495092, 0.10766892858911586, 0.09150339905675102, 0.10754153910848108, -0.0038021415630793837],
                             chord=1, color='b', N1=.5, N2=1, tol=0.0014553076272791395)

s = np.linspace(0, g_fit.arclength(1)[0], 100)
p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[1])
b_fit = beam_chen(g_fit, p, l, s, ignore_ends=True)
b_fit.g.internal_variables(b_fit.length)
b_fit.g.calculate_x1(b_fit.s)
b_fit.x = b_fit.g.x1_grid
b_fit.y = b_fit.g.x3(b_fit.x)


b_sol = beam_chen(g_sol, p, l, s, ignore_ends=True)
b_sol.g.internal_variables(b_sol.length)
b_sol.g.calculate_x1(b_fit.s)
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
plt.plot(b_sol.x, b_sol.y, 'k--', label='Sol')
plt.plot(b_fit.x, b_fit.y, 'k', label='Fit')
plt.legend()
plt.show()

plt.figure()
plt.plot(x, d_fit, 'b', label='dFit')
plt.plot(x, dd_fit, 'r', label='ddFit')
plt.plot(x, d_sol, 'b--', label='dSol')
plt.plot(x, dd_sol, 'r--', label='ddSol')
plt.ylim([-2, 2])
plt.legend()
plt.show()
