from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle

from aeropy.geometry.airfoil import CST, rotate
from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def cst(x, *A):
    b.g.D = A
    # b.g.internal_variables(b.g.length)
    y = []
    for xi in x:
        y.append(b.g.x3(np.array([xi]))[0])
    return y


plt.figure()

abaqus = np.loadtxt('case_study_6b_lower.csv', delimiter=',')
abaqus = abaqus.T
abaqus = abaqus[abaqus[:, 0].argsort()]
x = abaqus[:, 0]
y = abaqus[:, 1]

psi_spars = [0.2]
n = 2
p = 3
i = n*p+1
g_p = CoordinateSystem.pCST(D=i*[0.],
                            chord=[psi_spars[0], .7, .1],
                            color=['b', 'r', 'g'], N1=[1., 1., 1.], N2=[1., 1., 1.],
                            offset=-.05, continuity='C2', free_end=True,
                            root_fixed=True)
g = CoordinateSystem.pCST(D=i*[0.],
                          chord=[psi_spars[0], .7, .1],
                          color=['b', 'r', 'g'], N1=[1., 1., 1.], N2=[1., 1., 1.],
                          offset=-.05, continuity='C2', free_end=True,
                          root_fixed=True)
_, _, n_u = g._check_input([])
g.offset_s = 0
s = g.calculate_s([51, 51, 51], 1)
p = properties()
l = loads()
b = beam_chen(g, p, l, s)

popt, pcov = curve_fit(cst, x[:-1], y[:-1], p0=np.zeros(n_u), method='lm', maxfev=10000)
# b.g.D = [0.00854872,  0.00392494, -0.00053255, -0.0150316,  -0.01014882,  0.00011948]
# popt, pcov = curve_fit(cst, x, y, p0=np.array(
#     [0.00854872,  0.00392494, -0.00053255, -0.0150316,  -0.01014882,  0.00011948]))
print('Solution: ', popt)
print('Error: ', np.sqrt(np.diag(pcov)))
print('lower 1', b.g.cst[0].D)
print('lower 2', b.g.cst[1].D)
print('lower 3', b.g.cst[2].D)

# b.length = b.g.arclength(b.g.chord)
# print('length', b.length, b.g.arclength(b.g.chord))
# print('s', b.s)
b.g.calculate_x1(b.g.s)

x_fit = x
y_fit = []
for i in range(len(x_fit)):
    if i == (len(x_fit) - 1):
        y_fit.append(b.g.x3(np.array([b.g.cst[-1].offset_x + b.g.cst[-1].chord]))[0])
    else:
        y_fit.append(b.g.x3(np.array([x_fit[i]]))[0])
plt.scatter(x, y, c='b', label='Raw')
plt.scatter(x_fit, y_fit, c='r', label='Fit')

rho_p = b.g_p.radius_curvature(b.g_p.x1_grid, output_only=True)
rho = b.g.radius_curvature(b.g.x1_grid, output_only=True)
plt.figure()
plt.plot(b.s, b.p.young*b.p.inertia*(rho - rho_p), 'b')
plt.show()
