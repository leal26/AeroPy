from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings
import pickle

from aeropy.geometry.airfoil import CST, rotate
from aeropy.structural.beam import beam_chen, coupled_beams
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def format_u(input, g=None, g_p=None):
    return list(input)


def format_input(input, gu=None, gu_p=None, gl=None, gl_p=None):
    _, _, n_u = g_upper._check_input([])
    n_u = n_u - np.count_nonzero(gl.dependent)
    Au = format_u(input[:n_u], gu, gu_p)
    Al = format_u(input[n_u:], gl, gl_p)
    return Au, Al


def cst_coupled(x, *A):
    print(A)
    Au, Al = format_input(A, gl=a.bl.g)
    a.update(Au, Al)
    # b.g.internal_variables(b.g.length)
    y = []
    for xi in x[:index]:
        y.append(a.bl.g.x3(np.array([xi]))[0])
    for xi in x[index:]:
        y.append(a.bu.g.x3(np.array([xi]))[0])
    return y


warnings.filterwarnings("ignore", category=RuntimeWarning)

abaqus = np.loadtxt('case_study_8_lower.csv', delimiter=',')
abaqus = abaqus.T
abaqus = abaqus[abaqus[:, 0].argsort()]
xl = abaqus[:, 0]
yl = abaqus[:, 1]

abaqus = np.loadtxt('case_study_8_upper.csv', delimiter=',')
abaqus = abaqus.T
abaqus = abaqus[abaqus[:, 0].argsort()]
xu = abaqus[:, 0]
yu = abaqus[:, 1]

x = list(xl) + list(xu)
y = list(yl) + list(yu)

index = len(xl)

psi_spars = [0.2, 0.6, 0.9]
chords = []
for i in range(len(psi_spars)):
    if i == 0:
        chords.append(psi_spars[i])
    else:
        chords.append(psi_spars[i] - psi_spars[i-1])
chords.append(1-psi_spars[-1])

m = len(psi_spars)

n = 2
p = 4
i = n*p+1

g_upper = CoordinateSystem.pCST(D=i*[0, ],
                                chord=chords, color=['b', 'r', 'g', 'm'],
                                N1=len(chords)*[1.], N2=len(chords)*[1.],
                                offset=.05, continuity='C2', free_end=True,
                                root_fixed=True)
g_lower = CoordinateSystem.pCST(D=i*[0, ],
                                chord=chords, color=['b', 'r', 'g', 'm'],
                                N1=len(chords)*[1.], N2=len(chords)*[1.],
                                offset=-.05, continuity='C2', free_end=True,
                                root_fixed=True,
                                dependent=[True, True, True, False])

g_upper.calculate_s(N=[11, 11, 11, 11])
g_lower.calculate_s(N=[11, 11, 11, 11])
p_upper = properties()
p_lower = properties()
l_upper = loads(concentrated_load=[[-100*np.sqrt(2)/2, -100*np.sqrt(2)/2]], load_s=[1])
l_lower = loads(concentrated_load=[[100*np.sqrt(2)/2, 100*np.sqrt(2)/2]], load_s=[1-0.1])


a = coupled_beams(g_upper, g_lower, p_upper, p_lower, l_upper, l_lower, None,
                  None, ignore_ends=True, spars_s=psi_spars)

a.calculate_x()

_, _, n_u = g_upper._check_input([])
_, _, n_l = g_lower._check_input([])

popt, pcov = curve_fit(cst_coupled, x, y, p0=np.zeros(n_u+n_l - np.count_nonzero(a.bl.g.dependent)))
Au, Al = format_input(popt, gl=a.bl.g)
a.update(Au, Al)
print('Solution: ', popt)
print('upper 1', a.bl.g.cst[0].D)
print('upper 2', a.bl.g.cst[1].D)
print('upper 3', a.bl.g.cst[2].D)
print('upper 4', a.bl.g.cst[3].D)
print('lower 1', a.bl.g.cst[0].D)
print('lower 2', a.bl.g.cst[1].D)
print('lower 3', a.bl.g.cst[2].D)
print('lower 4', a.bl.g.cst[3].D)


# Plotting
plt.figure()
x_fit = xu
y_fit = []
for xi in x_fit:
    y_fit.append(a.bu.g.x3(np.array([xi]))[0])
plt.plot(xu, yu, 'b', label='Raw')
plt.plot(x_fit, y_fit, 'r--', label='Fit')

x_fit = xl
y_fit = []
for xi in x_fit:
    y_fit.append(a.bl.g.x3(np.array([xi]))[0])
plt.plot(xl, yl, 'b', label='Raw')
plt.plot(x_fit, y_fit, 'r--', label='Fit')

# y_fit = []
# for xi in x_fit:
#     y_fit.append(bl.g.x3(np.array([xi]), 'x1')[0])
# plt.plot(x_fit, y_fit, 'g--')
# plt.scatter([a.bl.g.spar_x], [a.bl.g.spar_y], c='g', label='Lower spar', zorder=20, s=40)

b = a.bu
rho_p = b.g_p.radius_curvature(b.g_p.x1_grid, output_only=True)
rho = b.g.radius_curvature(b.g.x1_grid, output_only=True)
plt.figure()
plt.plot(b.s, b.p.young*b.p.inertia*(rho - rho_p), 'b', label='upper')
b = a.bl
rho_p = b.g_p.radius_curvature(b.g_p.x1_grid, output_only=True)
rho = b.g.radius_curvature(b.g.x1_grid, output_only=True)
plt.plot(b.s, b.p.young*b.p.inertia*(rho - rho_p), 'r', label='lower')
plt.legend()
plt.show()
