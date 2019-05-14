import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

from aeropy.CST_2D.fitting import fitting_shape_coefficients, shape_parameter_study
from aeropy.xfoil_module import output_reader
from aeropy.geometry.airfoil import CST

directory = './'
filename = 'naca641212_upper.txt'
data_upper = np.genfromtxt(directory + filename)

n = 8
Data = shape_parameter_study(filename, n=n, solver='gradient', deltaz=0,
                             objective='squared_mean', surface='lower')
plt.figure()
x = np.linspace(2, 2*n, n)
plt.plot(x, Data['error'])
plt.scatter(x, Data['error'])
plt.grid()
plt.xlabel('Number of shape functions')
plt.ylabel('Error')
plt.show()
