import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem
from aeropy.CST_2D import calculate_c_baseline
from aeropy.geometry.airfoil import create_x, rotate

import cProfile
import pstats


def format_input(input, g=None, g_p=None):
    # def residual(deltaz):
    #     g.D[-1] = deltaz[0]
    #     g.internal_variables(target_length)
    #     rho_p = (1/g_p.chord)*dd_p/(1+d_p**2)**(3/2)
    #     dd = (2*n*g.D[-3] - 2*(g.N1+n)*g.D[-2])
    #     d = -g.D[-2] + g.D[-1]
    #     rho = (1/g.chord)*dd/(1+d**2)**(3/2)
    #     g.radius_curvature(np.array([g.chord]))
    #     rho_new = g_p.rho[0]
    #     print(deltaz, rho, rho_new, rho_p, rho_p_new)
    #     return(abs(rho_new-rho_p_new))

    # print('D before', g.D)
    g.D[:-1] = input
    error = 9999
    n = g.n - 2
    dd_p = (2*n*g_p.D[-3] - 2*(g_p.N1+n)*g_p.D[-2])
    d_p = -g_p.D[-2] + g_p.D[-1]
    rho_p = (1/g_p.chord)*dd_p/(1+d_p**2)**(3/2)
    target_length = g_p.arclength(g_p.chord)[0]
    # rho_p = dd_p/(1+d_p**2)**(3/2)
    # g.D[-1] = fsolve(residual, g.D[-1])[0]
    # BREAK
    # BREAK

    print('target', target_length)

    while error > 1e-4:
        before = g.D[-1]
        C = rho_p*g.chord
        dd = (2*n*g.D[-3] - 2*(g.N1+n)*g.D[-2])
        d = -g.D[-2] + g.D[-1]

        if C > 0:
            deltaz = np.sqrt((dd/C)**(2/3)-1) + g.D[-2]
        elif C < 0:
            deltaz = -np.sqrt((dd/C)**(2/3)-1) + g.D[-2]
        if not np.isnan(deltaz):
            g.D[-1] = deltaz
        else:
            break

        d = -g.D[-2] + g.D[-1]
        rho = (1/g.chord)*dd/(1+d**2)**(3/2)
        print('D after', g.D)
        print('terms', np.sqrt((dd/C)**(2/3)-1), g.D[-2])
        print('C', C, rho_p, g.chord, g_p.chord)
        print('dd', d, d_p, dd, dd_p)
        print('before', before)
        print('rho', rho, rho_p)

        g.internal_variables(target_length)
        after = g.D[-1]
        print('after', after)
        error = abs(after - before)
        print('deltaxi', g.D[-1], 'error', error)
        # BREAK
    return g.D


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

# rotated_abaqus = rotate({'x': abaqus_x, 'y': abaqus_y})

# plt.figure()
# plt.plot(abaqus_x, abaqus_y, 'b')
# plt.plot(rotated_abaqus['x'], rotated_abaqus['y'], 'r')
# plt.show()

# NACA0008
# higher order [0.1194, 0.0976, 0.1231, 0.0719, 0.1061, 0.1089, 0]
# N 100: 0.06430238065609895
# N 10: 0.011485823192467997
g = CoordinateSystem.CST(D=[0.1127, 0.1043, 0.0886, 0.1050, 0], chord=1,
                         color='b', N1=.5, N2=1, tol=0.0014553076272791395)
g_p = CoordinateSystem.CST(D=[0.1127, 0.1043, 0.0886, 0.1050, 0], chord=1,
                           color='k', N1=.5, N2=1, tol=0.0014553076272791395)

# g.x1_grid = create_x(1, n=100, distribution='polar')[::-1]

s = np.linspace(0, g.arclength(1)[0], 100)
# print('s', s)
# g.calculate_x1(s)
# g.darc = np.zeros(len(g.x1_grid))
# for i in range(len(g.x1_grid)):
#     s[i] = g.arclength_index(i)
#     print(i, g.arclength_index(i), g.arclength(g.x1_grid[i])[0])
#
# BREAK
p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[s[-1]], follower=True)

b = beam_chen(g, p, l, s, ignore_ends=True)
# b.g.D = [0.11621803608468839, 0.11164077707202186,
#          0.08581012060178156, 0.11474758149621056, -0.005855939703725668]
# b.g.D = [0.1126956897257928, 0.10811631223364697, 0.09120680993172335, 0.10763228117522525, 0]
# D = format_input([0.11621803608468839, 0.11164077707202186,
#                   0.08581012060178156, 0.11474758149621056], g=b.g, g_p=b.g_p)
# b.g.D = [0.11621803608468839, 0.11164077707202186,
#          0.08581012060178156, 0.11474758149621056, -0.005855939703725668]
# print('D', D)
# b.D = D
b.g_p.x1_grid = np.linspace(0, 1, 10)
print('Accurate', b.length)
print('Improper', b.g_p.improper_arclength_chord())
#
# b.g_p.internal_variables(b.length)
# b.g_p.calculate_x1(b.s)
# print('x', b.g_p.x1_grid)
BREAK


def test_function():
    b.g.internal_variables(b.length)
    b.g.calculate_x1(b.s)
    b.x = b.g.x1_grid
    b.y = b.g.x3(b.x)


# profile = cProfile.Profile()
# profile.runcall(test_function)
# ps = pstats.Stats(profile)
# ps.sort_stats('cumtime')
# ps.print_stats()
# BREAK
# b.parameterized_solver(format_input=format_input, x0=g.D[:-1])
b.g.internal_variables(b.length)
b.g.calculate_x1(b.s)
b.x = b.g.x1_grid
b.y = b.g.x3(b.x)


print('residual: ', b._residual(b.g.D))
# print('x', b.x)
# print('y', b.y)
# print('sizes', len(b.x), len(b.s))
# print('chord', b.g.chord)
g_p.calculate_x1(b.s)

plt.figure()
plt.plot(b.x, b.M/b.p.young/b.p.inertia, 'b')
plt.plot(b.x, b.g.rho, 'r')
plt.plot(b.x, b.g_p.rho, 'g')

plt.figure()
g_p.plot(label='Parent')
print('chord', b.g.chord)
plt.plot(b.x, b.y, '.5', label='Child: %.3f N' % -l.concentrated_load[0][-1], lw=3)
plt.scatter(abaqus_x, abaqus_y, c='.5', label='FEA: %.3f N' % -
            l.concentrated_load[0][-1], edgecolors='k', zorder=10, marker="^")

plt.legend()
plt.show()
