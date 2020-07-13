import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem
from aeropy.CST_2D import calculate_c_baseline
from aeropy.geometry.airfoil import create_x, rotate


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

# NACA0012
# higher order [0.1194, 0.0976, 0.1231, 0.0719, 0.1061, 0.1089, 0]
g = CoordinateSystem.CST(D=[0.1127, 0.1043, 0.0886, 0.1050, 0], chord=1,
                         color='b', N1=.5, N2=1, tol=0.0111414160997469)
g_p = CoordinateSystem.CST(D=[0.1127, 0.1043, 0.0886, 0.1050, 0], chord=1,
                           color='k', N1=.5, N2=1, tol=0.0111414160997469)

s = np.linspace(0, g.arclength(1)[0], 10)

p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[s[-1]], follower=True)

b = beam_chen(g, p, l, s, ignore_ends=True)
# b.g.D = [0.11621804, 0.11164078, 0.08581012, 0.11474758, -0.00585594]
# b.g.calculate_x1(b.s)
# b.y = b.g.x3(b.g.x1_grid)
# for i in range(len(b.s)):
#     print('%f, %f' % (b.g.x1_grid[i], b.y[i]))
# data = np.array([b.g.x1_grid, b.y]).T
# np.savetxt('points.csv', data, delimiter=",")

b.parameterized_solver(format_input=format_input, x0=g.D[:-1])
b.g.internal_variables(b.length)
b.g.calculate_x1(b.s)
b.x = b.g.x1_grid
b.y = b.g.x3(b.x)
print('residual: ', b._residual(b.g.D))
print('x', b.x)
print('y', b.y)
print('length', b.length, b.g.arclength(b.g.chord)[0])
g_p.calculate_x1(b.s)
g_p.plot(label='Parent')
plt.plot(b.x, b.y, '.5', label='Child: %.3f N' % -l.concentrated_load[0][-1], lw=3)
plt.scatter(rotated_abaqus['x'], rotated_abaqus['y'], c='.5', label='FEA: %.3f N' % -
            l.concentrated_load[0][-1], edgecolors='k', zorder=10, marker="^")
# plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.show()
