import matplotlib.pyplot as plt
import numpy as np
import json

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem

def format_input(input):
    return np.array(list(input[2:4]) + [0,0])

# Abaqus result for load [0, -1]
# abaqus_x = [-1.9931523e-39, 0.019805111, 0.039600205, 0.059363969, 0.079074532, 0.098731205, 0.11834421, 0.13791738, 0.15745261, 0.17695515, 0.19643463, 0.21590081, 0.23536223, 0.25482732, 0.27430427, 0.29380038, 0.31332213, 0.33287457, 0.35246155, 0.37208539, 0.39174658, 0.41144425, 0.4311758, 0.45093691, 0.47072202, 0.4905242, 0.51033556, 0.53014755, 0.54995131, 0.56973785, 0.58949906, 0.60922712, 0.62891579, 0.64856035, 0.66815758, 0.68770653, 0.70720816, 0.72666544, 0.74608386, 0.76547045, 0.78483433, 0.80418658, 0.82353973, 0.84290802, 0.86230612, 0.88174498, 0.901232, 0.92077732, 0.94040447, 0.96012044, 0.97989261, 0.99968362]
# abaqus_y = [-5.6866093e-41, 0.0005642015, 0.0014020122, 0.0027871644, 0.0047945748, 0.0072746095, 0.01008267, 0.013156435, 0.016460776, 0.019951705, 0.023572199, 0.02726537, 0.030982889, 0.03468034, 0.038314555, 0.041844603, 0.045229983, 0.048432998, 0.051417515, 0.05414886, 0.05659632, 0.058729947, 0.060523856, 0.061955255, 0.063003175, 0.06365294, 0.063891411, 0.063711293, 0.063111104, 0.062091313, 0.06066037, 0.058829125, 0.056613438, 0.054035675, 0.051119566, 0.047895364, 0.044395782, 0.040657114, 0.036721848, 0.032633409, 0.028439464, 0.02418977, 0.019940561, 0.015762938, 0.011728834, 0.0078962492, 0.004313848, 0.0010692433, -0.0016330082, -0.0035785499, -0.004843954, -0.0057649203]

abaqus_x = [-2.4561814e-36, 0.042830434, 0.081619181, 0.11706093, 0.15135017, 0.18578932, 0.22090997, 0.25701118, 0.29424605, 0.33458638, 0.37631872, 0.4194892, 0.46415922, 0.51042736, 0.55844355, 0.60836291, 0.65985292, 0.70932084, 0.74462074, 0.76366568, 0.77207059, 0.77436781, 0.77298594, 0.76925296, 0.76396841, 0.75766265, 0.75072503, 0.74347734, 0.73579228, 0.72851181, 0.72219682, 0.71769124, 0.7170409, 0.7248221, 0.7466259, 0.022113848, 0.062998489, 0.099586852, 0.13424191, 0.1685102, 0.20324115, 0.23882626, 0.27548063, 0.31424513, 0.35527575, 0.39772093, 0.44163153, 0.48708475, 0.5342049, 0.58315802, 0.63398248, 0.6854037, 0.72905624, 0.75586993, 0.76885015, 0.77379298, 0.77403677, 0.77135527, 0.76676732, 0.76091665, 0.75425178, 0.74712104, 0.73961163, 0.73207271, 0.72518462, 0.71964973, 0.7167412, 0.71869367, 0.73448193]
abaqus_y = [2.7112095e-36, 0.013644938, 0.044001155, 0.078382984, 0.11392248, 0.14931659, 0.18403351, 0.21772821, 0.25016478, 0.28260598, 0.31323448, 0.34179837, 0.36795166, 0.39115247, 0.41045329, 0.42400545, 0.42787951, 0.4141773, 0.38002855, 0.33443719, 0.28560939, 0.23607944, 0.1865052, 0.13704836, 0.087731533, 0.038534064, -0.010578668, -0.059646729, -0.11168896, -0.16378947, -0.21601416, -0.26842281, -0.32100561, -0.37273252, -0.40621749, 0.0055217547, 0.027570475, 0.060954697, 0.096128643, 0.13167301, 0.16678636, 0.20102447, 0.23411642, 0.26659808, 0.29816106, 0.32779297, 0.35520408, 0.37996736, 0.40137377, 0.41812, 0.42754599, 0.42381513, 0.39924994, 0.3579492, 0.31018928, 0.2608681, 0.21128196, 0.16175881, 0.11237314, 0.063119814, 0.013969492, -0.035115667, -0.085664503, -0.13772725, -0.18988574, -0.24218005, -0.29468915, -0.34748906, -0.38998169]
a4 = 3
a5 = 3
a2 = a4 + 2*a5
a3 = -2*a4-3*a5
g = CoordinateSystem.polynomial(D=[0, 0, a2, a3, a4, a5], chord = 1, color = 'b', n =6)
g_p = CoordinateSystem.polynomial(D=[0, 0, a2, a3, a4, a5], chord = 1, color = 'k', n =6)

s_tip = g.arclength(g.chord)[0]
s = np.linspace(0, s_tip, 41)

p = properties()
l = loads(concentrated_load = [[0, -100]], load_s = [s[-1]])
# s = np.linspace(0, 1, 10)
b = beam_chen(g, p, l, s)
# b.g.D = [0, 0, 3.03040918, -3.29131145, 0.55964773, -0.31282177]
b.g.calculate_x1(b.s)
b.y = b.g.x3(b.g.x1_grid)

np.savetxt("points.csv", np.array([b.g.x1_grid, b.y]).T, delimiter=",")

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

# b.g.D = [0, 0, 3, -3, 0, 0]
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