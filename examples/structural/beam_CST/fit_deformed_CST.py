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
    b.g.D = A + [-A[0]]
    b.g.internal_variables(b.g.length)
    y = b.g.x3(x)
    return y


plt.figure()

chord = 1

abaqus = np.loadtxt('lower_beam_hinge.csv', delimiter=',')
abaqus = abaqus.T
abaqus = abaqus[abaqus[:, 0].argsort()]
x = abaqus[:, 0]
y = abaqus[:, 1]
g_p = CoordinateSystem.CST(D=[0, 0, 0, 0, 0, 0],
                           chord=chord, color='k', N1=1, N2=1, deltaz=0, tol=1e-5, offset=-0.05)
g = CoordinateSystem.CST(D=[0, 0, 0, 0, 0, 0],
                         chord=chord, color='k', N1=1, N2=1, deltaz=0, tol=1e-5, offset=-0.05)
g.offset_s = 0
s = g.calculate_s(201, 1)
p = properties()
l = loads()
b = beam_chen(g, p, l, s)

popt, pcov = curve_fit(cst, x, y, p0=g.D[0:-1])
print('Solution: ', popt)
print('Error: ', np.sqrt(np.diag(pcov)))

b.length = b.g.arclength(b.g.chord)
print('length', b.length, b.g.arclength(b.g.chord))
print('s', b.s)
b.g.calculate_x1(b.s)

x_fit = np.linspace(0, b.g.chord, 100)
y_fit = b.g.x3(x_fit)
print('x', x_fit)
print('y', y_fit)
d = b.g.x3(np.array([b.g.chord]), diff='x1')
print('d', d)
print('c', b.g.chord)
print(b.g.deltaz)
plt.plot(x, y, 'b', label='Raw')
plt.plot(x_fit, y_fit, 'r--', label='Fit')

rho_p = b.g_p.radius_curvature(b.g_p.x1_grid, output_only=True)
rho = b.g.radius_curvature(b.g.x1_grid, output_only=True)
plt.figure()
plt.plot(b.s, b.p.young*b.p.inertia*(rho - rho_p), 'b')
plt.show()
