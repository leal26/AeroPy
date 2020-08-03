import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.optimize import fsolve

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem
from aeropy.CST_2D import calculate_c_baseline, dxi_u
from aeropy.geometry.airfoil import create_x, rotate, CST

# def con(input):
#     g.D = input
#     g.internal_variables(target_length)
#
#     dd = (2*n*g.D[-3] - 2*(g.N1+n)*g.D[-2])
#     d = -g.D[-2] + g.D[-1]
#     rho = (1./g.chord)*dd/(1+d**2)**(3/2)
#     print('rho', rho, rho_p, g.D[-1])
#     return(abs(rho-rho_p))


def con(input):
    g.D = format_input(input, g_p=g_p)
    g.internal_variables(target_length)

    lhs = -2/(g.D[0]**2*g.chord) + 2/(g_p.D[0]**2*g_p.chord)
    rhs = -1*g.chord/p.young/p.inertia
    return(lhs-rhs)


def format_input(input, g=None, g_p=None):
    g.deltaz = abaqus_y[-1]
    # n = len(g_p.D) - 2
    # Pn = g_p.D[-2]
    # Pn1 = g_p.D[-3]
    # Cn1 = input[-3]
    # Cn = input[-2]
    # dd_c = 2*n*Cn1-(1+2*n)*Cn
    # dd_p = 2*n*Pn1-(1+2*n)*Pn
    g.chord = abaqus_x[-1]

    d = np.zeros([2, 1])
    d[0] = g_p.x3(np.array([epsilon]))[0] - epsilon*b.g.deltaz
    d[1] = g_p.x3(np.array([epsilon]), diff='x1')[0] - b.g.deltaz/b.g.chord
    g.deltaz = abaqus_y[-1]/abaqus_x[-1]
    A0 = np.zeros(len(g.D)-1)
    D = np.zeros([2, 2])
    for i in range(2):
        A = np.copy(A0)
        A[i] = 1
        D[0, i] = CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=g.chord)
        D[1, i] = dxi_u(epsilon, A, 0, N1=0.5, N2=1)

    for i in range(2, 6):
        A = np.copy(A0)
        A[i] = input[i-2]
        d[0] -= CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=g.chord)
        d[1] -= dxi_u(epsilon, A, 0, N1=0.5, N2=1)
    A01 = np.linalg.solve(D, d)
    # print('inverse', A01)
    return [A01[0][0], A01[1][0]] + list(input) + [g.deltaz]


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

# rotated_abaqus = rotate({'x': abaqus_x, 'y': abaqus_y}, normalize=False)
# abaqus_x = rotated_abaqus['x']
# abaqus_y = rotated_abaqus['y']

epsilon = 0.01
g = CoordinateSystem.CST(D=[0.11397826, 0.10433884, 0.10241407, 0.10070566, 0.0836374, 0.11353368, 0.], chord=1,
                         color='b', N1=.5, N2=1)
g.name = 'proper integral'
g_p = CoordinateSystem.CST(D=[0.11397826, 0.10433884, 0.10241407, 0.10070566, 0.0836374, 0.11353368, 0.], chord=1,
                           color='k', N1=.5, N2=1)
g_p.name = 'proper integral'
s_epsilon = g_p.arclength(np.array([epsilon]))[0]
s = np.linspace(s_epsilon, g.arclength(np.array([1]))[0], 101)
# s = g_p.calculate_s(101, density='curvature')
p = properties()
l = loads(concentrated_load=[[0, -1]], load_s=[s[-1]])
print('s', s)
b = beam_chen(g, p, l, s, origin=epsilon)

n = g_p.n - 2
dd_p = (2*n*g_p.D[-3] - 2*(g_p.N1+n)*g_p.D[-2])
d_p = -g_p.D[-2] + g_p.D[-1]
rho_p = (1/g_p.chord)*dd_p/(1+d_p**2)**(3/2)
target_length = s[-1]

cons = {'type': 'eq', 'fun': con}
b.parameterized_solver(format_input=format_input, x0=g.D[2:-1])
# b.g.internal_variables(b.length)
b.g.calculate_x1(b.s, origin=b.origin, length_rigid=b.s[0])
b.x = b.g.x1_grid
b.y = b.g.x3(b.x)
b.g.radius_curvature(b.g.x1_grid)
b.g_p.radius_curvature(b.g_p.x1_grid)
# rotated_beam = rotate({'x': b.x, 'y': b.y}, normalize=False)
plt.figure()
plt.plot(b.x, b.M/b.p.young/b.p.inertia, c='b')
plt.plot(b.x, b.g.rho, c='r')
plt.plot(b.x, b.g_p.rho, c='g')
plt.figure()
plt.plot(b.x, b.M/b.p.young/b.p.inertia, c='b')
plt.plot(b.x, b.g.rho - b.g_p.rho, c='r')
plt.show()
print('residual: ', b._residual(b.g.D))
print(b.r)
print('x', b.x)
print('y', b.y)
print('length', target_length, b.length, b.g.arclength(b.g.chord)[0])
# g_p.calculate_x1(b.s)
b.g_p.plot(label='Parent')
plt.plot(b.x, b.y, '.5',
         label='Child: %.3f N' % -l.concentrated_load[0][-1], lw=3)
plt.scatter(abaqus_x, abaqus_y, c='.5', label='FEA: %.3f N' % -
            l.concentrated_load[0][-1], edgecolors='k', zorder=10, marker="^")
# plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.show()
