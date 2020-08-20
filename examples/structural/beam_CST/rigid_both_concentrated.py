import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.optimize import fsolve, fixed_point

from aeropy.structural.beam import beam_chen, airfoil
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem
from aeropy.CST_2D import calculate_c_baseline, dxi_u
from aeropy.geometry.airfoil import create_x, rotate, CST


def format_input(input, gu=None, gu_p=None, gl=None, gl_p=None):
    Au = format_input_upper(input[:-1], g=gu, g_p=gu_p)
    Al = format_input_lower(input[-1:], gu=gu, gu_p=gu_p, gl=gl, gl_p=gl_p)

    return Au, Al


def format_input_upper(input, g=None, g_p=None):
    def free_end(An_in):
        den = (1+(-An_in + g.zetaT)**2)**(1.5)
        An_out = (2*n*Cn1 - den*(g.chord/g_p.chord)*dd_p)/(1+2*n)
        return An_out
    error = 999
    while error > 1e-8:
        chord0 = np.copy(g.chord)

        # A5
        n = len(g_p.D) - 2
        Pn = g_p.D[-2]
        Pn1 = g_p.D[-3]
        Cn1 = input[-1]
        dd_p = (2*n*Pn1-(1+2*n)*Pn)/(1+(-Pn + g_p.zetaT)**2)**(1.5)
        An = float(fixed_point(free_end, g_p.D[-2]))
        temp = list(input) + [An]

        A_template = np.zeros(len(g.D)-1)
        # D
        D = np.zeros([2, 2])
        A = np.copy(A_template)
        A[0] = 1
        D[0, 0] = CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=g.chord)
        D[1, 0] = dxi_u(epsilon/g.chord, A, 0, N1=0.5, N2=1)
        D[0, 1] = epsilon
        D[1, 1] = 1

        # d
        d = np.zeros([2, 1])
        d[0] = g_p.x3(np.array([epsilon]))[0]
        d[1] = g_p.x3(np.array([epsilon]), diff='x1')[0]
        for i in range(1, n+1):
            A = np.copy(A_template)
            A[i] = temp[i-1]
            d[0] -= CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=g.chord)
            d[1] -= dxi_u(epsilon/g.chord, A, 0, N1=0.5, N2=1)

        det = D[0, 0]*D[1, 1] - D[0, 1]*D[1, 0]
        A0 = (1/det)*(D[1, 1]*d[0][0] - D[0, 1]*d[1][0])
        g.zetaT = (1/det)*(-D[1, 0]*d[0][0] + D[0, 0]*d[1][0])
        g.D = [A0] + list(temp) + [g.zetaT]
        g.internal_variables(a.bu.length, origin=epsilon)
        error = abs(g.chord-chord0)
    return [A0] + list(temp) + [g.zetaT]


def format_input_lower(A3, gu=None, gu_p=None, gl=None, gl_p=None):
    def free_end(An_in):
        den = (1+(-An_in + gl.zetaT)**2)**(1.5)
        An_out = (2*n*Cn1 - den*(gl.chord/gl_p.chord)*dd_p)/(1+2*n)
        dd_c1 = (1/gl.chord)*(2*n*Cn1-(1+2*n)*An_out)/(1+(-An_in + gl_p.zetaT)**2)**(1.5)
        dd_c2 = (1/gl.chord)*(2*n*Cn1-(1+2*n)*An_out)/(1+(-An_out + gl_p.zetaT)**2)**(1.5)
        print(dd_p, dd_c1, dd_c2)
        return An_out

    gl.D[3] = A3[0]
    gl.chord = gu.chord
    gl.zetaT = gu.zetaT
    gl.deltaz = gu.chord*gu.zetaT
    n = len(gl_p.D) - 2

    Pn = gl_p.D[-2]
    Pn1 = gl_p.D[-3]
    dd_p = (2*n*Pn1-(1+2*n)*Pn)/(1+(-Pn + gl_p.zetaT)**2)**(1.5)
    Cn1 = gl.D[-3]
    A0 = np.sqrt(gl_p.chord/gl.chord)*gl_p.D[0]
    An = float(fixed_point(free_end, gl.D[-2]))
    temp = [A3, An]

    A_template = np.zeros(len(gl.D)-1)
    # D
    D = np.zeros([2, 2])
    A = np.copy(A_template)
    A[1] = 1
    D[0, 0] = CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=gu.chord)
    D[1, 0] = dxi_u(epsilon/gu.chord, A, 0, N1=0.5, N2=1)
    A = np.copy(A_template)
    A[2] = 1
    D[0, 1] = CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=gu.chord)
    D[1, 1] = dxi_u(epsilon/gu.chord, A, 0, N1=0.5, N2=1)

    # d
    d = np.zeros([2, 1])
    d[0] = gl_p.x3(np.array([epsilon]))[0] - epsilon*gu.zetaT
    d[1] = gl_p.x3(np.array([epsilon]), diff='x1')[0] - gu.zetaT

    indexes = [0] + [*range(3, n+1)]
    for i in indexes:
        A = np.copy(A_template)
        if i == 0:
            A[0] = A0
        else:
            A[i] = temp[i-3]
        d[0] -= CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=gu.chord)
        d[1] -= dxi_u(epsilon/gu.chord, A, 0, N1=0.5, N2=1)

    det = D[0, 0]*D[1, 1] - D[0, 1]*D[1, 0]
    A1 = (1/det)*(D[1, 1]*d[0][0] - D[0, 1]*d[1][0])
    A2 = (1/det)*(-D[1, 0]*d[0][0] + D[0, 0]*d[1][0])
    gl.D[0] = A0
    gl.D[1] = A1
    gl.D[2] = A2
    gl.D[-2] = An
    gl.D = [A0, A1, A2, A3, An, gu.zetaT]

    return [A0, A1, A2, A3, An, gu.zetaT]


epsilon = 0.1
upper = np.loadtxt('upper_airfoil.csv', delimiter=',')
lower = np.loadtxt('lower_airfoil.csv', delimiter=',')
abaqus_x = list(upper[0, :]) + list(lower[0, :])
abaqus_y = list(upper[1, :]) + list(lower[1, :])
g_upper = CoordinateSystem.CST(D=[0.11415405, 0.10066617, 0.10787586, 0.08047336, 0.11218462, 0], chord=1,
                               color='b', N1=.5, N2=1)
g_lower = CoordinateSystem.CST(D=[-0.11415405, -0.10066617, -0.10787586, -0.08047336, -0.11218462, 0], chord=1,
                               color='b', N1=.5, N2=1)
g_upper.name = 'proper integral'
g_lower.name = 'proper integral'
g_p = CoordinateSystem.CST(D=[-0.11415405, -0.10066617, -0.10787586, -0.08047336, -0.11218462, 0], chord=1,
                           color='k', N1=.5, N2=1)
g_p.name = 'proper integral'
s_epsilon = g_p.arclength(np.array([epsilon]))[0]
s_upper = np.linspace(s_epsilon, g_p.arclength(np.array([1]))[0], 51)
s_lower = np.linspace(s_epsilon, g_p.arclength(np.array([1]))[0], 51)
p_upper = properties()
p_lower = properties()
l_upper = loads(concentrated_load=[[0, -1]], load_s=[s_upper[-1]])
l_lower = loads(concentrated_load=[[0, -1]], load_s=[s_lower[-1]])
a = airfoil(g_upper, g_lower, p_upper, p_lower, l_upper, l_lower, s_upper,
            s_lower, origin=epsilon, ignore_ends=True)
a.calculate_x()


# plt.show()
# BREAK
# print('s', a.bu.s)
# print('x', a.bu.g.x1_grid)
# target_length = s[-1]
#
a.parameterized_solver(format_input=format_input, x0=list(g_upper.D[1:-2]) + [g_lower.D[3]])
a.bl.g.radius_curvature(a.bl.g.x1_grid)
a.bu.g.radius_curvature(a.bu.g.x1_grid)
# b.gu.calculate_x1(b.s, origin=b.origin, length_rigid=b.s[0])
# a.bu.x = b.g.x1_grid
# b.y = b.g.x3(b.x)

plt.figure()
plt.plot(a.bu.g.x1_grid, a.bu.g.rho - a.bu.g_p.rho, label='Upper')
plt.plot(a.bl.g.x1_grid, a.bl.g.rho - a.bl.g_p.rho, label='Lower')
print('rho child', a.bl.g.rho)
print('rho parent', a.bl.g_p.rho)
plt.legend()

plt.figure()
plt.plot(a.bu.g_p.x1_grid, a.bu.g_p.x3(a.bu.g_p.x1_grid), 'b',
         label='Upper Parent', lw=3)
plt.plot(a.bu.g.x1_grid, a.bu.g.x3(a.bu.g.x1_grid), '.5',
         label='Upper Child: %.3f N' % -l_upper.concentrated_load[0][-1], lw=3)
plt.plot(a.bl.g_p.x1_grid, a.bl.g_p.x3(a.bl.g_p.x1_grid), 'b', linestyle='dashed',
         label='Lower Parent', lw=3)
plt.plot(a.bl.g.x1_grid, a.bl.g.x3(a.bl.g.x1_grid), '.5', linestyle='dashed',
         label='Lower Child: %.3f N' % -l_upper.concentrated_load[0][-1], lw=3)

plt.scatter(abaqus_x, abaqus_y, c='.5', label='FEA: %.3f N' % -
            l_upper.concentrated_load[0][-1], edgecolors='k', zorder=10, marker="^")
plt.legend()
plt.show()
