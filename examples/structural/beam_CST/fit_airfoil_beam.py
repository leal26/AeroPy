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


chord = 1

abaqus = np.loadtxt('case_study_C1b.csv', delimiter=',')
abaqus = abaqus.T
abaqus = abaqus[abaqus[:, 0].argsort()]
x = abaqus[:, 0]
y = abaqus[:, 1]
psi_spars = [0.2]

D_n4 = [0.288197, 0.37918722, 0.21903394, 0.25925831, 0.24259769, 0.27278782,
        0.06528661, 0.07283469, 0.06692653, 0.0559131, -0.01268226, -0.00134043,
        -0.08743061, -0.0013288]

D_n3 = [3.06132404e-01, 3.31838944e-01, 2.44074856e-01, 2.30082260e-01,
        2.71512524e-01, 4.98187521e-02, 8.59903075e-02, 4.89862447e-02,
        1.22835553e-04, -7.05141464e-02, -1.95302773e-02]

D_n2 = [0.3208729, 0.28457915, 0.2431048, 0.27140013, 0.06470603, 0.06222951,
        -0.02810757, -0.04412038]

Db_n3 = [0.30645304, 0.33092568, 0.2454664,  0.23047635, 0.27141484, 0.04881688,
         0.09722178, 0.02404006]

# g_p = CoordinateSystem.pCST(D=Db_n3,
#                             chord=[.2, .8], color=['b', 'r'],
#                             N1=[.5, 1], N2=[1, 1], continuity='C2',
#                             free_end=True, rigid_LE=True)
# g = CoordinateSystem.pCST(D=Db_n3,
#                           chord=[.2, .8], color=['b', 'r'],
#                           N1=[.5, 1], N2=[1, 1], continuity='C2',
#                           free_end=True, rigid_LE=True)

g_p = CoordinateSystem.pCST(D=D_n2,
                            chord=[.2, .7, .1], color=['b', 'r', 'g'],
                            N1=[.5, 1, 1], N2=[1, 1, 1], continuity='C2',
                            free_end=True, rigid_LE=True)
g = CoordinateSystem.pCST(D=D_n2,
                          chord=[.2, .7, .1], color=['b', 'r', 'g'],
                          N1=[.5, 1, 1], N2=[1, 1, 1], continuity='C2',
                          free_end=True, rigid_LE=True)
_, _, n_u = g._check_input([])
g.offset_s = 0
# s = g.calculate_s([51, 51], 1)
s = g.calculate_s([51, 51, 51], 1)
p = properties()
l = loads()
b = beam_chen(g, p, l, s)

b.g.calculate_x1(b.g.s)
x_fit = b.g.x1_grid
dy = []
ddy = []

for xi in x_fit:
    dy.append(b.g.x3(np.array([xi]), 'x1')[0])
    ddy.append(b.g.x3(np.array([xi]), 'x11')[0])
plt.figure()
plt.plot(x_fit, dy, 'b--', label='dy (parent)')
plt.plot(x_fit, ddy, 'r--', label='ddy (parent)')
plt.legend()

print(b.g.cst[0].D)
print(b.g.cst[1].D)
popt, pcov = curve_fit(cst, x[:], y[:], p0=[0.0*i for i in range(n_u)], maxfev=10000)
# b.g.D = [0.10688066710180918, 0.04881687846018838, 0.0972217791492235]
# print('last')
# b.g.D = popt
print('Solution: ', popt)
print('Error: ', np.sqrt(np.diag(pcov)))
print('D1', b.g.cst[0].D)
print('D2', b.g.cst[1].D)
# print('D3', b.g.cst[2].D)
# b.length = b.g.arclength(b.g.chord)
# print('length', b.length, b.g.arclength(b.g.chord))
# print('s', b.s)
b.g.calculate_x1(b.g.s)

x_fit = b.g.x1_grid
dy = []
ddy = []

for xi in x_fit:
    dy.append(b.g.x3(np.array([xi]), 'x1')[0])
    ddy.append(b.g.x3(np.array([xi]), 'x11')[0])
plt.plot(x_fit, dy, 'b', label='dy (child)')
plt.plot(x_fit, ddy, 'r', label='ddy (child)')
plt.legend()

x_fit = x
y_fit = []
for xi in x_fit:
    y_fit.append(b.g.x3(np.array([xi]))[0])
plt.figure()
plt.scatter(x, y, c='b', label='Raw')
plt.scatter(x_fit, y_fit, c='r', label='Fit')
plt.legend()

rho_p = b.g_p.radius_curvature(b.g_p.x1_grid, output_only=True)[:]
rho = b.g.radius_curvature(b.g.x1_grid, output_only=True)[:]
plt.figure()
plt.plot(b.g.x1_grid[:], b.p.young*b.p.inertia*(rho - rho_p), 'b')


plt.show()
