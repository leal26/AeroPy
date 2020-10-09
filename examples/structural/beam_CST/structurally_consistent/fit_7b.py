from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings
import pickle

from aeropy.geometry.airfoil import CST, rotate
from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def cst_upper(x, *A):
    bu.g.D = A
    # b.g.internal_variables(b.g.length)
    y = []
    for xi in x:
        y.append(bu.g.x3(np.array([xi]))[0])
    return y


def cst_lower(x, *A):
    # print('lower', [A0, A1, A2, A3, A4, A5, A6, A7, A8])
    bl.g.D = A

    bl.g.cst[0].chord = bl.g.spar_x[0]
    current_length = bl.g.cst[0].arclength(bl.g.cst[0].chord)
    target_length = bl.g.cst[0].length
    penalization = abs(target_length-current_length)
    print(penalization, target_length, current_length)

    y = []
    for xi in x:
        y.append(bl.g.x3(np.array([xi]))[0])
    return y


warnings.filterwarnings("ignore", category=RuntimeWarning)

plt.figure()

chord = 1
abaqus = np.loadtxt('case_study_7b_upper.csv', delimiter=',')
abaqus = abaqus.T
abaqus = abaqus[abaqus[:, 0].argsort()]
xu = abaqus[:, 0]
yu = abaqus[:, 1]

psi_spars = [0.2]
gu_p = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0., 0., 0.],
                             chord=[psi_spars[0], 1-psi_spars[0]],
                             color=['b', 'r'], N1=[1., 1.], N2=[1., 1.],
                             offset=.05, continuity='C2', free_end=True,
                             root_fixed=True)
gu = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0., 0., 0.],
                           chord=[psi_spars[0], 1-psi_spars[0]],
                           color=['b', 'r'], N1=[1., 1.], N2=[1., 1.],
                           offset=.05, continuity='C2', free_end=True,
                           root_fixed=True)
_, _, n_u = gu._check_input([])
gu.offset_s = 0
su = gu.calculate_s([51, 51], 1)
p = properties()
l = loads()
bu = beam_chen(gu, p, l, su)


popt, pcov = curve_fit(cst_upper, xu, yu, p0=np.zeros(n_u))
bu.g.D = popt
print('Solution upper: ', popt)
# bu.g.D = [0.0085628,  0.01159214, 0.01473726, 0.0176863,  0.06739102, 0.0548096]

bu.g.calculate_x1(bu.g.s)


abaqus = np.loadtxt('case_study_7b_lower.csv', delimiter=',')
abaqus = abaqus.T
abaqus = abaqus[abaqus[:, 0].argsort()]
xl = abaqus[:, 0]
yl = abaqus[:, 1]
gl_p = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                             chord=[psi_spars[0], 0.7, 0.1],
                             color=['b', 'r', 'g'], N1=[1., 1., 1.], N2=[1., 1., 1.],
                             offset=-.05, continuity='C2', free_end=True,
                             root_fixed=False, dependent=[True, False, False])
gl = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                           chord=[psi_spars[0], 0.7, 0.1],
                           color=['b', 'r', 'g'], N1=[1., 1., 1.], N2=[1., 1., 1.],
                           offset=-.05, continuity='C2', free_end=True,
                           root_fixed=True, dependent=[True, False, False])
gl.g_independent = bu.g
_, _, n_l = gl._check_input([])
print('n', n_l)
gl.offset_s = 0
sl = gl.calculate_s([51, 51, 51], 1)
p = properties()
l = loads()
bl = beam_chen(gl, p, l, su)


popt, pcov = curve_fit(cst_lower, xl[:-1], yl[:-1], p0=np.zeros(n_l), maxfev=10000)
bl.g.D = popt
print('Solution lower: ', popt)
# bl.g.D = [5.25846692e-03,  2.72634735e-03, -3.53410506e-04, -1.70242793e-02,
#           -1.34587434e-02, -1.01267466e-02, -1.38736304e-05,  2.96782298e-05]
bl.g.calculate_x1(bl.g.s)
# bl.g._calculate_spar(0, debug=True)
print('lower 1', bl.g.cst[0].D)
print('lower 2', bl.g.cst[1].D)
print('lower 3', bl.g.cst[2].D)
print('spars', bl.g.spar_x, bl.g.spar_y)
x_fit = xu
y_fit = []
for xi in x_fit:
    y_fit.append(bu.g.x3(np.array([xi]))[0])
plt.plot(xu, yu, 'b', label='Raw')
plt.plot(x_fit, y_fit, 'r--', label='Fit')

x_fit = xl
y_fit = []
i = 0
for xi in x_fit:
    if i == (len(x_fit) - 1):
        y_fit.append(bl.g.x3(np.array([bl.g.cst[-1].offset_x + bl.g.cst[-1].chord]))[0])
    else:
        y_fit.append(bl.g.x3(np.array([xi]))[0])
    i += 1
plt.plot(xl, yl, 'b', label='Raw')
plt.plot(x_fit, y_fit, 'r--', label='Fit')
# for i in range(len(x_fit)):
#     print(x_fit[i], y_fit[i])
print('total chord', bl.g.cst[-1].offset_x + bl.g.cst[-1].chord)
x_fit = xu
y_fit = []
for xi in x_fit:
    y_fit.append(bu.g.x1(np.array([xi]), 'theta1')[0])
plt.plot(x_fit, y_fit, 'b-.', label='sy')

print('total chord', bl.g.cst[-1].offset_x + bl.g.cst[-1].chord)
x_fit = xu
y_fit = []
for xi in x_fit:
    y_fit.append(-bu.g.x3(np.array([xi]), 'theta1')[0])
plt.plot(x_fit, y_fit, 'r-,', label='sx')

upper_spar_x = bl.g.spar_psi_upper*bu.g.cst[0].chord
upper_spar_y = bl.g.spar_xi_upper*bu.g.cst[0].chord
FEA_spar_x = [0.199995, 0.199984]
FEA_spar_y = [-0.0517113, 0.0482887]
plt.scatter([bl.g.spar_x, upper_spar_x], [bl.g.spar_y, upper_spar_y],
            c='g', label='Lower spar', zorder=20, s=40)
plt.scatter(FEA_spar_x, FEA_spar_y, c='b', s=40)
print('x/psi outside', upper_spar_x, bl.g.spar_psi_upper)
print('FEA distance', np.sqrt((FEA_spar_x[0]-FEA_spar_x[1])**2 +
                              (FEA_spar_y[0]-FEA_spar_y[1])**2))
print('CST distance', np.sqrt((upper_spar_x-bl.g.spar_x)**2 +
                              (upper_spar_y-bl.g.spar_y)**2))
b = bu
rho_p = b.g_p.radius_curvature(b.g_p.x1_grid, output_only=True)
rho = b.g.radius_curvature(b.g.x1_grid, output_only=True)
plt.figure()
plt.plot(b.s, b.p.young*b.p.inertia*(rho - rho_p), 'b', label='upper')
b = bl
rho_p = b.g_p.radius_curvature(b.g_p.x1_grid, output_only=True)
rho = b.g.radius_curvature(b.g.x1_grid, output_only=True)
plt.plot(b.s, b.p.young*b.p.inertia*(rho - rho_p), 'r', label='lower')
plt.legend()
plt.show()
