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

    g = CoordinateSystem.CST(D=list(A) + [-A[0]], chord=1., N1=1., N2=1., color='b')
    print(g.D)
    y = []
    for xi in x:
        y.append(g.x3(np.array([xi]))[0])
    return y


n = 4

plt.figure()

raw_airfoil = np.loadtxt('input_A7.csv', delimiter=',')

raw_airfoil = raw_airfoil[raw_airfoil[:, 0].argsort()]
x = raw_airfoil[:, 0]
y = raw_airfoil[:, 1]

popt, pcov = curve_fit(cst, x, y, p0=(n+1)*[0], maxfev=10000)
print('Solution: ', popt)
print('Error: ', np.sqrt(np.diag(pcov)))


g = CoordinateSystem.CST(D=list(popt) + [-popt[0]], chord=1., N1=1., N2=1., color='b')
s = g.calculate_s(N=51)
p = properties()
l = loads()
b = beam_chen(g, p, l, s)
print('D', b.g.D)
b.g.calculate_x1(b.g.s)

x_fit = x
y_fit = []
ddy_fit = []
for xi in x_fit:
    y_fit.append(b.g.x3(np.array([xi]))[0])
    ddy_fit.append(b.g.x3(np.array([xi]), diff='x11')[0])
plt.scatter(x, y, c='b', label='Raw')
plt.scatter(x_fit, y_fit, c='r', label='Fit')
plt.plot(x_fit[1:], ddy_fit[1:], 'g', label='dd')
plt.legend()
plt.show()
