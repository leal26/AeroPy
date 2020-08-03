from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
from scipy.optimize import fsolve

from aeropy.geometry.airfoil import CST, rotate
from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem
from aeropy.CST_2D import dxi_u


def calculate_A(input):
    d = np.zeros([2, 1])

    d[0] = b.g_p.x3(np.array([epsilon]))[0] - epsilon*b.g.deltaz
    d[1] = b.g_p.x3(np.array([epsilon]), diff='x1')[0] - b.g.deltaz/b.g.chord

    A0 = np.zeros(len(g.D)-1)
    D = np.zeros([2, 2])
    for i in range(2):
        A = np.copy(A0)
        A[i] = 1
        D[0, i] = CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=b.g.chord)
        D[1, i] = dxi_u(epsilon, A, 0, N1=0.5, N2=1)

    for i in range(2, 6):
        A = np.copy(A0)
        A[i] = input[i-2]
        d[0] -= CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=b.g.chord)
        d[1] -= dxi_u(epsilon, A, 0, N1=0.5, N2=1)
    A01 = np.linalg.solve(D, d)
    print(A01)
    return [A01[0][0], A01[1][0]]


def calculate_chord(input):
    A4 = input[-2]
    A5 = input[-1]
    n = 5
    Pn = b.g_p.D[-2]
    Pn1 = b.g_p.D[-3]
    Cn1 = A4
    Cn = A5
    dd_c = 2*n*Cn1-(1+2*n)*Cn
    dd_p = 2*n*Pn1-(1+2*n)*Pn
    return b.g_p.chord*dd_c/dd_p


def f(deltaz):
    b.g.deltaz = deltaz[0]
    s_current = b.g.arclength(b.g.chord, origin=epsilon)[0]
    print(deltaz, b.length, s_current)
    return b.length - s_current


def cst(x, A2, A3, A4, A5, deltaz):

    input = [A2, A3, A4, A5]

    b.g.chord = calculate_chord(input)
    b.g.deltaz = deltaz*b.g.chord
    b.g.D = calculate_A(input) + list(input) + [b.g.deltaz]
    # b.g.deltaz = fsolve(f, b.g.deltaz, xtol=1e-4)[0]
    # b.g.D[-1] = b.g.deltaz
    y = b.g.x3(x)
    return y


# def cst(x, A0, A1, A2, A3, A4, A5):
#     b.g.D = [A0, A1, A2, A3, A4, A5, -0.0058551501/0.99970293]
#     b.deltaz = -0.0058551501/0.99970293
#     # b.length = b.g.arclength(chord=chord)[0]
#     # b.g.internal_variables(b.length)
#     y = b.g.x3(x)
#     return y


# Airfoil (NACA0012) under 1 N
abaqus_x = [0, 0.0046626413, 0.018402023, 0.048702486, 0.087927498, 0.12745513,
            0.16709827, 0.20679522, 0.24651785, 0.28625035, 0.3259826,
            0.36570814, 0.40542319, 0.44512612, 0.48481697, 0.52449685,
            0.5641672, 0.60382956, 0.64348471, 0.68313259, 0.72277194,
            0.76240051, 0.80201453, 0.84160942, 0.88117975, 0.92072034,
            0.96022779, 0.99970293]
abaqus_y = [0, 0.0080771167, 0.015637144, 0.024132574, 0.030406451,
            0.034420185, 0.037082944, 0.038775206, 0.039692041, 0.039953224,
            0.039647505, 0.038851682, 0.037638135, 0.036076538, 0.034232546,
            0.032165118, 0.029923303, 0.027543247, 0.025045833, 0.022435322,
            0.019699335, 0.016810443, 0.013729585, 0.010411576, 0.0068127746,
            0.0029010926, -0.0013316774, -0.0058551501]

x = abaqus_x
y = abaqus_y


g = CoordinateSystem.CST(D=[0.11397826, 0.10433884, 0.10241407, 0.10070566, 0.0836374,  0.11353368, 0],
                         chord=1, color='k', N1=.5, N2=1, deltaz=0, tol=None)
g.name = 'proper integral'
epsilon = 0.02

s_epsilon = g.arclength(np.array([epsilon]))[0]
s = np.linspace(s_epsilon, g.arclength(np.array([1]))[0], 101)

p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[s[-1]])

b = beam_chen(g, p, l, s, origin=epsilon)

popt, pcov = curve_fit(cst, x, y, p0=b.g_p.D[2:])
print('Solution: ', popt)
print('Error: ', np.sqrt(np.diag(pcov)))

b.g.chord = calculate_chord(list(popt)[:-1])
b.g.deltaz = popt[-1]*b.g.chord
b.g.D = calculate_A(popt) + list(popt)

# b.g.D[-1] = float(fsolve(f, b.g.deltaz, xtol=1e-4)[0])
# b.g.deltaz = b.g.D[-1]

print('deltaz', b.g.deltaz)
b.g.calculate_x1(b.s, origin=b.origin, length_rigid=b.s[0])
b.g_p.calculate_x1(b.s, origin=b.origin, length_rigid=b.s[0])
print(b.g.x1_grid)
print('Parent', b.g_p.chord, max(b.g_p.x1_grid))
print('child', b.g.chord, max(b.g.x1_grid))
print('arclength', b.g.arclength(b.g_p.chord, origin=epsilon)[0], b.g.arclength(
    b.g.chord, origin=epsilon)[0], b.g.arclength(max(b.g.x1_grid), origin=epsilon)[0])
x_fit = b.g.x1_grid
y_fit = b.g.x3(b.g.x1_grid)
print('D', b.g.D)
# b._residual(b.g.D)

b.g.radius_curvature(b.g.x1_grid)
print('Second derivative')
n = 5
dd_c = 2*n*b.g.D[-3]-(1+2*n)*b.g.D[-2]
dd_p = 2*n*b.g_p.D[-3]-(1+2*n)*b.g_p.D[-2]
print(dd_p, b.g_p.x3(b.g_p.x1_grid[-1], diff='x11'))
print(dd_c, b.g.x3(b.g.x1_grid[-1], diff='x11'))

plt.figure()
plt.plot(x, y, 'b', label='Raw')
plt.plot(x_fit, y_fit, 'r--', label='Fit')

plt.figure()
plt.plot(b.g.x1_grid, b.g.rho, 'b', label='Child')
plt.plot(b.g_p.x1_grid, b.g_p.rho, 'r', label='Parent')

plt.figure()
plt.plot(b.g.x1_grid, b.M/b.p.young/b.p.inertia, 'b', label='Child')
plt.plot(b.g.x1_grid, b.g.rho - b.g_p.rho, 'r', label='Parent')
plt.show()
