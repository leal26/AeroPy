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


g_5 = CoordinateSystem.CST(D=[0.11397826, 0.10433884, 0.10241407, 0.10070566, 0.0836374,  0.11353368, 0], chord=1,
                           color='b', N1=.5, N2=1, deltaz=0)
g = CoordinateSystem.CST(D=[0.1127, 0.1043, 0.0886, 0.1050, 0], chord=1,
                         color='b', N1=.5, N2=1, deltaz=0)

x = np.linspace(0, 1, 1000)
plt.figure()
for i in range(6):
    d = dxi_u(x, np.ones(6), 0, N1=0.5, N2=1, indexes=[i])
    plt.plot(x, d, label=i)
plt.legend()
plt.ylabel('Influence of Shape component on $dy/dx$')

plt.figure()
for i in range(6):
    D_i = [0]*6
    D_i[i] = 1
    d = ddxi_u(x, D_i, N1=0.5, N2=1)
    plt.plot(x, d, label=i)
plt.legend()
plt.ylabel('Influence of Shape component on $ddy/ddx$')

plt.figure()
for i in range(6):
    D_i = [0]*6
    D_i[i] = 1
    d = CST(x, Au=D_i, deltasz=0, N1=0.5, N2=1, c=1)
    plt.plot(x, d, label=i)
plt.legend()
plt.ylabel('Influence of Shape component on $y$')
plt.show()
D_fit = [0.12023167906378877, 0.10112629883623848, 0.12518045636490516,
         0.07661460343152796, 0.10804192536758454, 0.11319177161706924, 0]
D_sol = [0.11398155727411789, 0.10708199288796764, 0.1048666057536085,
         0.10382292685834452, 0.08645752995277559, 0.11617334120517299, 0]

s = g.calculate_s(101, density='curvature')
# s = np.linspace(0, g.arclength(np.array([1]))[0], 101)
# s = g.calculate_s(201, density='curvature')
p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[1])

b_fit = beam_chen(g_5, p, l, s)
# b_fit.g.D = D_fit
b_fit.g.internal_variables(b_fit.length)
b_fit.g.calculate_x1(b_fit.s)
b_fit.x = b_fit.g.x1_grid
b_fit.y = b_fit.g.x3(b_fit.x)


b_sol = beam_chen(g_5, p, l, s)
b_sol.g.D = D_sol
b_sol.g.internal_variables(b_sol.length)
b_sol.g.calculate_x1(b_fit.s)
b_sol.x = b_sol.g.x1_grid
b_sol.y = b_sol.g.x3(b_sol.x)

x = b_fit.g.x1_grid
d_fit = b_fit.g.x3(x, diff='x1')
dd_fit = b_fit.g.x3(x, diff='x11')
d_sol = b_sol.g.x3(x, diff='x1')
dd_sol = b_sol.g.x3(x, diff='x11')

x_ds_fit = b_fit.g.x1(x, diff='theta1')
x_dds_fit = b_fit.g.x1(x, diff='theta11')
y_ds_fit = b_sol.g.x3(x, diff='theta1')
y_dds_fit = b_sol.g.x3(x, diff='theta11')

plt.figure()
# plt.plot(x, x_ds_fit, 'b', label='x_ds')
# plt.plot(x, x_dds_fit, 'r', label='x_dds')
# plt.plot(x, y_ds_fit, 'b--', label='y_ds')
# plt.plot(x, y_dds_fit, 'r--', label='y_dds')
plt.plot(x, d_fit, 'b', label='d_fit')
plt.plot(x, dd_fit, 'r', label='dd_fit')
plt.plot(x, d_sol, 'b--', label='d_sol')
plt.plot(x, dd_sol, 'r--', label='dd_sol')
plt.ylim([-2, 2])
plt.legend()
# plt.show()

rho_p_o = np.copy(b_sol.g_p.radius_curvature(b_sol.g_p.x1_grid, output_only=True))
rho_p_p = np.copy(b_sol.g_p.radius_curvature(b_sol.g_p.x1_grid, output_only=True, parametric=True))

rho_o = np.copy(b_sol.g.radius_curvature(b_sol.g.x1_grid, output_only=True))
rho_p = np.copy(b_sol.g.radius_curvature(b_sol.g.x1_grid, output_only=True, parametric=True))
plt.figure()
plt.plot(rho_o - rho_p_o, label='Original')
plt.plot(rho_p - rho_p_p, label='Parametric')
# plt.plot(, 'r')
plt.legend()
plt.show()
BREAK

print('x', x)
print('dfit', d_fit)
print(5*D_sol[-3])
plt.figure()
plt.plot(b_sol.x, b_sol.y, 'k--', label='Sol')
plt.plot(b_fit.x, b_fit.y, 'k', label='Fit')
plt.legend()

plt.figure()
plt.plot(x, d_fit, 'b', label='dFit')
plt.plot(x, dd_fit, 'r', label='ddFit')
plt.plot(x, d_sol, 'b--', label='dSol')
plt.plot(x, dd_sol, 'r--', label='ddSol')
plt.ylim([-2, 2])
plt.legend()

b_sol.g.radius_curvature(b_sol.g.x1_grid)
b_fit.g.radius_curvature(b_fit.g.x1_grid)
plt.figure()
plt.plot(b_sol.x, b_sol.g.rho, 'b', label='Sol')
plt.scatter(b_fit.x, b_fit.g.rho, c='g', label='Fit')
plt.ylabel('Curvature')
plt.legend()

plt.figure()
plt.plot(b_sol.x, b_sol.g.rho - b_sol.g_p.rho, 'b', label='Sol')
plt.plot(b_fit.x, b_fit.g.rho - b_fit.g_p.rho, 'g--', label='Fit')
plt.ylabel('LHS')
plt.legend()
plt.show()
