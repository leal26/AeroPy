from aeropy.geometry.parametric import CoordinateSystem
from aeropy.structural.stable_solution import (properties, boundary_conditions, euler_bernoulle)
from aeropy.structural.shell import shell
from aeropy.xfoil_module import output_reader

from scipy.integrate import quad
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from optimization_tools.DOE import DOE
import pickle

def input_function(x):
    return np.array([0,0] + list(x))

def find_chord(P):
    def _to_minimize(l):
        def _to_integrate(y):
            den = np.sqrt(1-P**2/E**2/I**2*(l*y-y**2/2)**2)
            if np.isnan(den):
                return 100
            else:
                return 1/den
        l = l[0]
        current_L = quad(_to_integrate, 0, l)[0]
        return abs(L-current_L)

    return minimize(_to_minimize, L, method = 'Nelder-Mead',).x[0]

def find_deflection(y, l, P):
    def _to_integrate(y):
        num = P/E/I*(l*y-y**2/2)
        den = np.sqrt(1-P**2/E**2/I**2*(l*y-y**2/2)**2)
        if np.isnan(den):
            return 100
        else:
            return num/den

    x = []
    for y_i in y:
        x_i = current_L = quad(_to_integrate, 0, y_i)[0]
        x.append(x_i)
    return x

def dsdy(y, l):
    den = np.sqrt(1-P**2/E**2/I**2*(l*y-y**2/2)**2)
    return 1/den

def dxdy(y, l):
    num = P/E/I*(l*y-y**2/2)
    den = np.sqrt(1-P**2/E**2/I**2*(l*y-y**2/2)**2)
    return num/den

np.set_printoptions(precision=5)

# Results from Abaqus
abaqus_data = pickle.load(open('neutral_line.p', 'rb'))

# Beam properties
bp = properties()
bc = boundary_conditions(concentrated_load=np.array([[0.0, 0.0, -1.0], ]))
EB_solution = bc.concentrated_load[0][2]/(6*bp.young*bp.inertia) * \
    np.array([0, 0, 3, -1])

curve_parent = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord = 1, color = 'b')
curve_child = CoordinateSystem.polynomial(D=EB_solution, chord = 1, color  ='0.5')

beam = shell(curve_parent, curve_child, bp, bc)
chord_bounds = [[0., 2],]
beam.g_p.bounds = chord_bounds
beam.g_c.bounds = chord_bounds
beam.theta1 = np.linspace(0, beam.g_p.arclength()[0], 10)
beam.update_parent()
beam.update_child()
eulerBernoulle = beam.g_c.r(beam.g_c.x1_grid)

# Find stable solution
bounds = np.array(((-0.02,0.02),
                  (-0.02,0.02)))


# beam.minimum_potential(x0=[0,0], input_function = lambda x: [0,0] + list(x),
                       # bounds = bounds)

# Euler Bernoulle
EB = euler_bernoulle(bp, -1, 'concentrated', beam.g_c, x=beam.g_c.x1_grid)
EB.analytical_solutions()
EB.bending_strain()
EB.free_energy()
EB.strain_energy()
print('U_EB: ', EB.U)
print('U_PL: ', beam.U)
# mod = 1/(1-0.3**2)
# print('U_corrected: ', beam.U/mod)
#
print('Bending_strain_EB: ', EB.B)
print('Bending_strain_PL: ', beam.rho)
#
# print('Phi_EB: ', EB.phi)
# print('Phi_PL: ', beam.phi)
# print('Phi_PL_B: ', beam.phi_B)
#
# print('ddx3_EB: ', EB.g.x3(beam.g_c.x1_grid, 'x11'))
# print('ddx3_PL: ', beam.g_c.x3(beam.g_c.x1_grid, 'x11'))

print('theta1',beam.g_c.x1(beam.g_c.x1_grid, 'theta1'))
print('theta11', beam.g_c.x1(beam.g_c.x1_grid, 'theta11'))
P = -1
E = bp.young
L = 1
I = bp.inertia
chord = find_chord(P)
print('chord_chen: ', chord)
print('chord_PL(c): ', beam.g_c.chord)
print('chord_PL(p): ', beam.g_p.chord)
# print('a')
# print(beam.g_c.a)
# print('da')
# print(beam.g_c.da)
# print('dA')
# print(beam.g_c.dA)
# beam_chen_x = np.linspace(0,chord,10)
# beam_chen_y = find_deflection(beam_chen_x, chord, P)
#
# [x,y,z] = eulerBernoulle.T
# beam.g_p.plot(label='Parent')
# plt.plot(x,z, 'k', lw = 3, label='Euler-Bernoulle', linestyle = '-', zorder=0)
# beam.g_c.plot(label='Child', linestyle = '--')
# plt.plot(beam_chen_x, beam_chen_y, 'r', label='Chen', linestyle = '--', lw = 3)
# plt.scatter(abaqus_data['coord'][0:401:40,0], abaqus_data['coord'][0:401:40,1], c='g', label='FEA', edgecolors='k', zorder = 10)
# plt.legend()
# plt.show()
#
#
#
# dxdy_chen = dxdy(beam_chen_x, chord)
# dsdy_chen = dsdy(beam_chen_x, chord)
# dxdy_leal = beam.g_c.x3(beam.g_c.x1_grid, 'x1')
# dsdy_leal = []
#
# for x1_i in beam.g_c.x1_grid:
#     dr = beam.g_c.r(np.array([x1_i]), 'x1')
#     dsdy_leal.append(np.sqrt(np.inner(dr, dr)[0,0]))
# plt.figure()
# plt.plot(beam_chen_x, dxdy_chen, 'r', linestyle = '-', label='Chen')
# plt.plot(beam.g_c.x1_grid, dxdy_leal, 'b', linestyle = '-', label='Leal')
# plt.title('dx3/dx1')
# plt.legend()
#
# plt.figure()
# plt.plot(beam_chen_x, dsdy_chen, 'r', linestyle = '-', label='Chen')
# plt.plot(beam.g_c.x1_grid, dsdy_leal, 'b', linestyle = '-', label='Leal')
# plt.title('ds/dx1')
# plt.legend()
# plt.show()
#
plt.figure()
plt.plot(beam.g_c.x1_grid,EB.B,'r')
plt.plot(beam.g_c.x1_grid,-beam.rho[0][0],'b')
plt.title('curvature radius')
plt.show()

plt.figure()
plt.plot(beam.g_c.x1_grid,EB.B,'r')
plt.plot(beam.g_c.x1_grid,beam.g_c.x3(beam.g_c.x1_grid, 'x11'),'b')
plt.title('ddx3/ddx1')
plt.show()
#
# plt.figure()
# plt.plot(beam.g_c.x1_grid,EB.phi,'r')
# plt.plot(beam.g_c.x1_grid,beam.width*beam.phi,'b')
# plt.title('Free energy')
# plt.show()
