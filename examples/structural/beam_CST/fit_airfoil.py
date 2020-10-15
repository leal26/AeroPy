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
    g = CoordinateSystem.pCST(D=A,
                              chord=[psi_spars[0], .7, .1],
                              color=['b', 'r', 'g'], N1=[.5, 1., 1.], N2=[1., 1., 1.],
                              offset=-.0, continuity='C2')
    y = []
    for xi in x:
        y.append(g.x3(np.array([xi]))[0])
    return y


plt.figure()

raw_airfoil = np.loadtxt('sc1095_upper.dat', delimiter=' ', skiprows=1)

raw_airfoil = raw_airfoil[raw_airfoil[:, 0].argsort()]
x = raw_airfoil[:, 0]
y = raw_airfoil[:, 1]

psi_spars = [0.2]

n = 4
p = 3
i = n*p+2
popt, pcov = curve_fit(cst, x, y, p0=np.zeros(i), maxfev=10000)
print('Solution: ', popt)
print('Error: ', np.sqrt(np.diag(pcov)))

g = CoordinateSystem.pCST(D=popt,
                          chord=[psi_spars[0], .7, .1],
                          color=['b', 'r', 'g'], N1=[.5, 1., 1.], N2=[1., 1., 1.],
                          offset=-.0, continuity='C2')
s = g.calculate_s([51, 51, 51], 1)
p = properties()
l = loads()
b = beam_chen(g, p, l, s)
print('D1', g.cst[0].D)
print('D2', g.cst[1].D)
print('D3', g.cst[2].D)
# b.length = b.g.arclength(b.g.chord)
# print('length', b.length, b.g.arclength(b.g.chord))
# print('s', b.s)
b.g.calculate_x1(b.g.s)

x_fit = b.g.x1_grid
y_fit = []
for xi in x_fit:
    y_fit.append(b.g.x3(np.array([xi]))[0])
plt.plot(x, y, 'b', label='Raw')
plt.plot(x_fit, y_fit, 'r--', label='Fit')
plt.show()
