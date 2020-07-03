from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle

from aeropy.geometry.airfoil import CST
from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def cst(x, A0, A1, A2, A3):
    b.g.D = [A0, A1, A2, A3, -A0]
    b.g.internal_variables(b.length)
    # b.g.calculate_x1(b.s)
    # b.g.chord = b.g.x1_grid[-1]
    # print('c1', b.g.chord)
    # b.g.chord = 0.889322917
    # b.g.deltaz = -A0*b.g.chord
    y = b.g.x3(x)
    # print('c2', b.g.chord)
    return y


plt.figure()
# EB
# abaqus_x = [0.,  0.11111111, 0.22222222, 0.33333333, 0.44444444,
#             0.55555556, 0.66666667, 0.77777778, 0.88888889, 1.]
# abaqus_y = [0., -0.0001019, -0.00039193, -0.00084656, -0.00144229, -
#             0.00215559, -0.00296296, -0.00384088, -0.00476582, -0.00571429]

# EB (distributed)
# abaqus_x = [0.,  0.11111111, 0.22222222, 0.33333333, 0.44444444,
#             0.55555556, 0.66666667, 0.77777778, 0.88888889, 1.]
# abaqus_y = [0.00000000e+00, -4.90996582e-05, -1.82028001e-04, -3.79188713e-04,
#             -6.23598319e-04, -9.00886189e-04, -1.19929453e-03, -1.50967840e-03,
#             -1.82550569e-03, -2.14285714e-03]
# abaqus_x = [-6.143523e-37, 0.049518712, 0.097414128, 0.14182274, 0.18338385, 0.22280164, 0.26078638, 0.29787692, 0.33455449, 0.3694379, 0.40472263, 0.44081238, 0.47817007, 0.51734918,
#             0.55898702, 0.60366577, 0.65133238, 0.70004117, 0.74728715, 0.78860354, 0.82396096, 0.85445607, 0.8812418, 0.90515256, 0.92684138, 0.94658399, 0.96465063, 0.98153377, 0.99782264]
# abaqus_y = [-1.0867065e-36, 0.0069836318, 0.025579287, 0.051580817, 0.08195062, 0.11507001, 0.14982793, 0.18554173, 0.22168092, 0.25601763, 0.2899411, 0.32300428, 0.35462314, 0.38394535,
#             0.40963244, 0.4294889, 0.44009688, 0.43760237, 0.42023277, 0.39134949, 0.35533783, 0.31509525, 0.27228072, 0.22779004, 0.18217251, 0.13567843, 0.088506527, 0.040896866, -0.0069233831]
# # Rao
abaqus_x = [0, 0.131510417, 0.2109375, 0.294270833, 0.384114583, 0.460286458,
            0.529947917, 0.598958333, 0.67578125, 0.746744792, 0.815104167, 0.889322917]
# y_B1 = [0, -0.008624495, -0.025178451, -0.043081887, -0.063603111, -0.081562353, -0.109830402, -0.143237635, -0.177948689, -0.221685023, -0.256406224, -0.298869177, -0.349063735, -0.390248233, -0.462349008]
abaqus_y = [0, -0.014878806, -0.03765024, -0.068739984, -0.108779046, -0.149511719, -
            0.193474559, -0.235286058, -0.281631611, -0.328230469, -0.375198317, -0.422476963]

# abaqus_data = pickle.load(open('neutral_line.p', 'rb'))
# abaqus_x = abaqus_data['coord'][0:401:40, 0]
# abaqus_y = abaqus_data['coord'][0:401:40, 1]
x = np.array(abaqus_x)
y = np.array(abaqus_y)

g = CoordinateSystem.CST(D=[0, 0, 0, 0, 0], chord=1, color='k', N1=1, N2=1, deltaz=0)

s = np.linspace(0, 1, 10)
p = properties()
l = loads()
b = beam_chen(g, p, l, s)
# BRAKE
# Fit
popt, pcov = curve_fit(cst, x, y, p0=[0, 0, 0, 0])
print('Solution: ', popt)
print('Error: ', np.sqrt(np.diag(pcov)))
b.g.D = np.append(popt, -popt[:1])
b.g.internal_variables(b.length)
print('c0', b.g.chord)
# b.g.chord = 0.889322917
# print('c1', b.g.chord)
# b.g.deltaz = b.g.D[-1]*b.g.chord
# print('deltaz1', b.g.deltaz, -0.422476963)
print('length', b.length, b.g.arclength(b.g.chord))
# b.g.calculate_x1(b.s)
print('s', b.s)
b.g.calculate_x1(b.s)
x_fit = b.g.x1_grid
y_fit = b.g.x3(x_fit)
print('x', x_fit)
print('y', y_fit)
print('c', b.g.chord)
print(b.g.deltaz)
plt.plot(x, y, 'b', label='Raw')
plt.plot(x_fit, y_fit, 'r--', label='Fit')
plt.show()
