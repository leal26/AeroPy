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

    # b.g.internal_variables(b.g.length)
    y = []
    for xi in x:
        y.append(bl.g.x3(np.array([xi]))[0])
    return y


warnings.filterwarnings("ignore", category=RuntimeWarning)

plt.figure()

chord = 1
abaqus = np.loadtxt('case_study_6b_upper.csv', delimiter=',')
abaqus = abaqus.T
abaqus = abaqus[abaqus[:, 0].argsort()]
xu = abaqus[:, 0]
yu = abaqus[:, 1]

psi_spars = [0.2]
gu_p = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0.],
                             chord=[psi_spars[0], 1-psi_spars[0]],
                             color=['b', 'r'], N1=[1., 1.], N2=[1., 1.],
                             offset=.05, continuity='C2', free_end=True,
                             root_fixed=True)
gu = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0.],
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
bu.g.calculate_x1(bu.g.s)


abaqus = np.loadtxt('case_study_6b_lower.csv', delimiter=',')
abaqus = abaqus.T
abaqus = abaqus[abaqus[:, 0].argsort()]
xl = abaqus[:, 0]
yl = abaqus[:, 1]
gl_p = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0., 0., 0.],
                             chord=[psi_spars[0], 0.7, 0.1],
                             color=['b', 'r', 'g'], N1=[1., 1., 1.], N2=[1., 1., 1.],
                             offset=-.05, continuity='C2', free_end=True,
                             root_fixed=False)
gl = CoordinateSystem.pCST(D=[0., 0., 0., 0., 0., 0., 0., 0.],
                           chord=[psi_spars[0], 0.7, 0.1],
                           color=['b', 'r', 'g'], N1=[1., 1., 1.], N2=[1., 1., 1.],
                           offset=-.05, continuity='C2', free_end=True,
                           root_fixed=True)
gl.g_independent = bu.g
_, _, n_l = gl._check_input([])
print('n', n_l)
gl.offset_s = 0
sl = gl.calculate_s([51, 51, 51], 1)
p = properties()
l = loads()
bl = beam_chen(gl, p, l, su)


popt, pcov = curve_fit(cst_lower, xl, yl, p0=np.zeros(n_l))
bl.g.D = popt
print('Solution lower: ', popt)
bl.g.calculate_x1(bl.g.s)

print('lower 1', bl.g.cst[0].D)
print('lower 2', bl.g.cst[1].D)
print('lower 3', bl.g.cst[2].D)
x_fit = xu
y_fit = []
for xi in x_fit:
    y_fit.append(bu.g.x3(np.array([xi]))[0])
plt.plot(xu, yu, 'b', label='Raw')
plt.plot(x_fit, y_fit, 'r--', label='Fit')

x_fit = xl
y_fit = []
for xi in x_fit:
    y_fit.append(bl.g.x3(np.array([xi]))[0])
plt.plot(xl, yl, 'b', label='Raw')
plt.plot(x_fit, y_fit, 'r--', label='Fit')
y_fit = []
for xi in x_fit:
    y_fit.append(bl.g.x3(np.array([xi]), 'x1')[0])

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
