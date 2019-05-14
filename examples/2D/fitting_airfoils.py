from scipy.optimize import minimize
from scipy.spatial.distance import directed_hausdorff
import matplotlib.pyplot as plt
import numpy as np
import math

from aeropy.CST_2D import CST
from aeropy.CST_2D.fitting import fitting_shape_coefficients, shape_parameter_study

directory = './'
filename = 'naca641212_upper.txt'
data_upper = np.genfromtxt(directory + filename)
filename = 'naca641212_lower.txt'
data_lower = np.genfromtxt(directory + filename)

deltaz, Al = fitting_shape_coefficients(data_lower, n=5, surface='lower',
                                        solver='gradient',
                                        objective='squared_mean')
deltaz, Au = fitting_shape_coefficients(data_upper, n=5, surface='upper',
                                        solver='gradient',
                                        objective='squared_mean')
print(Au)
print(Al)
y_u = CST(data_upper[:, 0], 1, deltasz=0, Au=Au)
y_l = CST(data_lower[:, 0], 1, deltasz=0, Al=Al)
plt.figure()
plt.scatter(data_upper[:, 0], data_upper[:, 1], label='raw_upper')
plt.scatter(data_lower[:, 0], data_lower[:, 1], label='raw_lower')
plt.plot(data_upper[:, 0], y_u, label='upper')
plt.plot(data_lower[:, 0], y_l, label='lower')
plt.legend()
plt.show()
