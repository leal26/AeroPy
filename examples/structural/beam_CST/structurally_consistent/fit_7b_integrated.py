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


def cst_all(x, *A):

    Au = A[:n_u]
    Al = A[n_u:]

    xu = x[:int(len(x)/2)]
    xl = x[int(len(x)/2):]

    yu = cst_upper(xu, Au)
    bl.g.g_independent = bu.g
    yl = cst_lower(xl, Al)

    # current_length = bu.g.cst[0].arclength(bu.g.cst[0].chord)
    # target_length = bu.g.cst[0].length
    # penalization_u = abs(target_length-current_length)
    #
    # bl.g.cst[0].chord = bl.g.spar_x[0]
    current_length = bl.g.cst[0].arclength(bl.g.cst[0].chord)
    target_length = bl.g.cst[0].length
    penalization_l = abs(target_length-current_length)
    # print(abs(target_length-current_length), penalization_u, penalization_l)
    # calculated_spar_x = bl.g.cst[0].chord
    # calculated_spar_y = bl.g.x3(np.array([calculated_spar_x]))
    # dx = abs(bl.g.spar_x - 0.199995)
    # dy = abs(bl.g.spar_y - -0.0517113)
    # d = 1000*np.sqrt(dx**2+dy**2)
    # print('x', bl.g.spar_x, calculated_spar_x, dx)
    # print('y', bl.g.spar_x, calculated_spar_y, dy)
    yu = np.array(yu)
    yl = np.array(yl) + penalization_l

    return(list(yu) + list(yl))


def cst_upper(x, *A):
    bu.g.D = A[0]
    # b.g.internal_variables(b.g.length)
    y = []
    for xi in x:
        y.append(bu.g.x3(np.array([xi]))[0])
    return y


def cst_lower(x, *A):
    # print('lower', [A0, A1, A2, A3, A4, A5, A6, A7, A8])
    bl.g.D = A[0]

    # b.g.internal_variables(b.g.length)
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

_, _, n_l = gl._check_input([])

gl.offset_s = 0
sl = gl.calculate_s([51, 51, 51], 1)
p = properties()
l = loads()
bl = beam_chen(gl, p, l, su)

popt, pcov = curve_fit(cst_all, list(xu[:-1]) + list(xl[:-1]), list(yu[:-1]) + list(yl[:-1]),
                       p0=np.zeros(n_u + n_l), maxfev=100000)
bu.g.D = popt[:n_u]
bl.g.g_independent = bu.g
bl.g.D = popt[n_u:]
print('Solution lower: ', popt)
bl.g.calculate_x1(bl.g.s)
bu.g.calculate_x1(bu.g.s)

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
for i in range(len(x_fit)):
    print(x_fit[i], y_fit[i])
print('total chord', bl.g.cst[-1].offset_x + bl.g.cst[-1].chord)
y_fit = []
for xi in x_fit:
    y_fit.append(bl.g.x3(np.array([xi]), 'x1')[0])

upper_spar_x = bu.g.cst[0].chord
upper_spar_y = bu.g.x3(np.array([upper_spar_x]))
plt.scatter([bl.g.spar_x, upper_spar_x], [bl.g.spar_y, upper_spar_y],
            c='g', label='Calculated spar', zorder=20, s=40)
plt.scatter([0.199995, 0.199984], [-0.0517113, 0.0482887], c='b', s=40)

b = bu
rho_p = b.g_p.radius_curvature(b.g_p.x1_grid, output_only=True)
rho = b.g.radius_curvature(b.g.x1_grid, output_only=True)
plt.figure()
plt.plot(b.g.x1_grid, b.p.young*b.p.inertia*(rho - rho_p), 'b', label='upper')
b = bl
rho_p = b.g_p.radius_curvature(b.g_p.x1_grid, output_only=True)
rho = b.g.radius_curvature(b.g.x1_grid, output_only=True)
plt.plot(b.g.x1_grid, b.p.young*b.p.inertia*(rho - rho_p), 'r', label='lower')
plt.legend()
plt.show()
