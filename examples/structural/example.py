import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem

def format_input(input):
    return np.array(list(input[2:4]) + [0,0])

# Abaqus result for load [0, -1]
abaqus_x = [-6.143523e-37, 0.049518712, 0.097414128, 0.14182274, 0.18338385, 0.22280164, 0.26078638, 0.29787692, 0.33455449, 0.3694379, 0.40472263, 0.44081238, 0.47817007, 0.51734918, 0.55898702, 0.60366577, 0.65133238, 0.70004117, 0.74728715, 0.78860354, 0.82396096, 0.85445607, 0.8812418, 0.90515256, 0.92684138, 0.94658399, 0.96465063, 0.98153377, 0.99782264]
abaqus_y = [-1.0867065e-36, 0.0069836318, 0.025579287, 0.051580817, 0.08195062, 0.11507001, 0.14982793, 0.18554173, 0.22168092, 0.25601763, 0.2899411, 0.32300428, 0.35462314, 0.38394535, 0.40963244, 0.4294889, 0.44009688, 0.43760237, 0.42023277, 0.39134949, 0.35533783, 0.31509525, 0.27228072, 0.22779004, 0.18217251, 0.13567843, 0.088506527, 0.040896866, -0.0069233831]

g = CoordinateSystem.polynomial(D=[0, 0, 3, -3, 0, 0], chord = 1, color = 'b', n =6)
g_p = CoordinateSystem.polynomial(D=[0, 0, 3, -3, 0, 0], chord = 1, color = 'k', n =6)

s_tip = g.arclength(g.chord)[0]
s = np.linspace(0, s_tip, 21)

p = properties()
l = loads(concentrated_load = [[0, -1]], load_s = [s[-1]])
# s = np.linspace(0, 1, 10)
b = beam_chen(g, p, l, s)


# b.g.D = [0, 0, 3.03040918, -3.29131145, 0.55964773, -0.31282177]
# b.g_p.calculate_x1(b.s)
# b.g.calculate_x1(b.s)
# b.g.radius_curvature(b.g.x1_grid)
#
# plt.subplot(4, 1, 1)
# plt.ylabel('y')
# y_p = b.g_p.x3(b.g_p.x1_grid)
# y = b.g.x3(b.g_p.x1_grid)
# plt.plot(b.g_p.x1_grid, y_p, 'k', label='y')
# plt.plot(b.g.x1_grid, y, '0.5',label='y')
#
# plt.subplot(4, 1, 2)
# plt.ylabel('dy')
# y_p = b.g_p.x3(b.g_p.x1_grid, diff='x1')
# y = b.g.x3(b.g.x1_grid, diff='x1')
# plt.plot(b.g_p.x1_grid, y_p, 'k', label='y')
# plt.plot(b.g.x1_grid, y, '0.5',label='y')
#
# plt.subplot(4, 1, 3)
# plt.ylabel('ddy')
# y = b.g_p.x3(b.g_p.x1_grid, diff='x11')
# y = b.g.x3(b.g.x1_grid, diff='x11')
# plt.plot(b.g_p.x1_grid, y_p, 'k', label='ddy')
# plt.plot(b.g.x1_grid, y, '0.5', label='ddy')
#
# plt.subplot(4, 1, 4)
# plt.ylabel('curvature')
# plt.plot(b.g_p.x1_grid, b.g_p.rho, 'k', label='curvature')
# plt.plot(b.g.x1_grid, b.g.rho, '0.5',label='curvature')
# plt.show()

# b.g.D = [0, 0, 3.03040918, -3.29131145, 0.55964773, -0.31282177]
# b.g.calculate_x1(b.s)
# b.y = b.g.x3(b.g.x1_grid)
# plt.figure()
# plt.plot(b.g.x1_grid, b.y, 'g', label='Fit')
# dy = b.g.x3(b.g.x1_grid, diff='x1')
# ddy = b.g.x3(b.g.x1_grid, diff='x11')
#
# x1= b.g.x1_grid
# D = b.g.D
# print('D', D)
# dy_check = 5*D[5]*x1**4 + 4*D[4]*x1**3 + 3*D[3]*x1**2 + 2*D[2]*x1 + D[1]
# ddy_check = 20*D[5]*x1**3 + 12*D[4]*x1**2 + 6*D[3]*x1 + 2*D[2]
# print('x', b.g.x1_grid)
# print('dy', dy)
# print('ddy', ddy)
#
# plt.plot(b.g.x1_grid, dy, 'g--', label='dFit')
# plt.plot(b.g.x1_grid, dy_check, 'k--', label='dCheck')
# plt.plot(b.g.x1_grid, ddy, 'g-.', label='ddFit')
# plt.plot(b.g.x1_grid, ddy_check, 'k-.', label='ddCheck')
# plt.show()

b.g.D = [0, 0, 3, -3, 0, 0]
# for i in range(len(b.x)):
#     print('%f, %f' % (b.x[i], b.y[i]))
b.parameterized_solver(format_input = format_input)

g_p.calculate_x1(s)
g_p.plot(label='Parent')
plt.plot(b.x, b.y, '.5', label='Child: %.3f N' % -l.concentrated_load[0][1], lw = 3)
plt.scatter(abaqus_x, abaqus_y)
# plt.ylim([-0.01, 0.04])
plt.legend()
# plt.gca().set_aspect('equal', adjustable='box')
plt.show()
