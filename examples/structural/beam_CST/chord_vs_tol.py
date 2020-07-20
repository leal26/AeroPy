import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem
from aeropy.CST_2D import calculate_c_baseline
from aeropy.geometry.airfoil import create_x, rotate

import cProfile
import pstats


def calculate_chord(tol, N):
    g.tol = tol
    g_p.tol = tol

    s = np.linspace(0, g.arclength(1)[0], N)
    l = loads()
    g.s = s

    b = beam_chen(g, p, l, s, ignore_ends=True)
    b.g.D = [0.1127, 0.1043, 0.0886, 0.1050, 0]
    b.g.internal_variables(b.length)
    b.g.calculate_x1(b.s)
    b.x = b.g.x1_grid
    # print(np.sort(b.g.x1_grid))
    return(max(b.x), b.g.chord, b.g.x1_grid == np.sort(b.g.x1_grid))


def residual(tol):
    max_x, chord, valid = calculate_chord(tol, N)
    return(abs(chord - max_x))


def format_input(input, g=None):
    # print('input', input)
    # nd_deltaz = input[-1]
    # # COnsidering BC for zero derivative at the root
    # error = 9999
    # c_P = 1
    # Au_P = [0.1127, 0.1043, 0.0886, 0.1050]
    # AC_u0 = Au_P[0]
    # deltaz_C = nd_deltaz*c_P
    # while error > 1e-9:
    #     before = AC_u0
    #     Au_C = [AC_u0] + list(input[:-1])
    #     # print('Au', Au_C)
    #     c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz_C)
    #     deltaz_C = nd_deltaz*c_C
    #     AC_u0 = np.sqrt(c_P/c_C)*Au_P[0]
    #     error = abs(AC_u0 - before)
    # return list(Au_C) + [nd_deltaz]
    return list(input) + [0]


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

rotated_abaqus = rotate({'x': abaqus_x, 'y': abaqus_y})

# plt.figure()
# plt.plot(abaqus_x, abaqus_y, 'b')
# plt.plot(rotated_abaqus['x'], rotated_abaqus['y'], 'r')
# plt.show()

# NACA0008
# higher order [0.1194, 0.0976, 0.1231, 0.0719, 0.1061, 0.1089, 0]
g = CoordinateSystem.CST(D=[0.1127, 0.1043, 0.0886, 0.1050, 0], chord=1,
                         color='b', N1=.5, N2=1)
g_p = CoordinateSystem.CST(D=[0.1127, 0.1043, 0.0886, 0.1050, 0], chord=1,
                           color='k', N1=.5, N2=1)

# g.x1_grid = create_x(1, n=100, distribution='polar')[::-1]


# g.calculate_x1(s)
# g.darc = np.zeros(len(g.x1_grid))
# for i in range(len(g.x1_grid)):
#     s[i] = g.arclength_index(i)
#     print(i, g.arclength_index(i), g.arclength(g.x1_grid[i])[0])
#
# BREAK
p = properties()


orders_N = np.linspace(1, 2, 2)
Ns = np.around(10**orders_N)
orders = np.linspace(-5, -1, 20)
tolerances = 10**orders

orders_tolerance = []
optimal_tolerance = []
maximum_tolerance = []
N_list = []
for j in range(len(Ns)):
    N = int(Ns[j])
    print(j, N)
    tol = fsolve(residual, 1e-3)
    max_x, chord, valid = calculate_chord(tol, N)
    if valid.all():
        orders_tolerance.append(N)
        optimal_tolerance.append(tol[0])
        N_list.append(N)
    print(N, tol[0], valid.all(), max_x, chord)

data = np.array([N_list, optimal_tolerance]).T
pickle.dump(data, open('tolerances.p', 'wb'))

plt.figure()
plt.plot(orders_tolerance, optimal_tolerance, 'k')
plt.xlabel('Order of magnitude')
plt.legend()
plt.show()


# N = 100
# # print(g.arclength(1)[0])
# chords = np.zeros(len(tolerances))
# max_x = np.zeros(len(tolerances))
# for i in range(len(tolerances)):
#
#     tol = tolerances[i]
#     max_xi, chord, valid = calculate_chord(tol, N)
#     print(tol, max_xi, g.arclength(1)[0], valid.all())
#     chords[i] = chord
#     max_x[i] = max_xi
# print('max_x', max_x)
# f = interp1d(max_x, tolerances)
# # print('Optimal tolerance: ', f(1))
# plt.figure()
# plt.plot(tolerances, chords, 'b', label='chord')
# plt.plot(tolerances, max_x, 'r', label='max_x')
# plt.xlabel('Tolerance')
# plt.legend()
# plt.show()
