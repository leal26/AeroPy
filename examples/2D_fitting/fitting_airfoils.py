from scipy.optimize import minimize
from scipy.spatial.distance import directed_hausdorff
import matplotlib.pyplot as plt
import numpy as np
import pickle
import math
import os

from aeropy.CST_2D import CST
from aeropy.CST_2D.fitting import fitting_shape_coefficients, shape_parameter_study
from aeropy.geometry.airfoil import rotate

f = {'names': [], 'upper': [], 'lower': [], 'Au': [], 'Al': [], 'du': [], 'dl': []}

logfile = './log.txt'
log = open(logfile, 'w')
directory = "./airfoils/"
file = "uag8814320.dat"

filename = os.path.join(directory, file)
extracted_data = False
n_skip = 0
while not extracted_data:
    try:
        data = np.genfromtxt(filename, skip_header=1, skip_footer=n_skip,
                             encoding="latin1")
        extracted_data = True
    except(ValueError):
        n_skip += 1
data = rotate({'x': data[:, 0], 'y': data[:, 1]},
              move_to_origin=True, both_surfaces=True)
data = np.array([data['x'], data['y']]).T
i_break = np.where(data[:, 0] == 0.)[0][0]
data_upper = data[0:i_break+1, :]
data_lower = data[i_break:, :]

deltaz_l, Al = fitting_shape_coefficients(data_lower, n=4, surface='lower',
                                          solver='differential_evolution',
                                          objective='squared_mean',
                                          deltaz=0)
deltaz_u, Au = fitting_shape_coefficients(data_upper, n=4, surface='upper',
                                          solver='differential_evolution',
                                          objective='squared_mean',
                                          deltaz=0)
f['names'].append(file[:-4])
f['upper'].append(data_upper)
f['lower'].append(data_lower)
f['Au'].append(Au)
f['Al'].append(Al)
f['du'].append(deltaz_u)
f['dl'].append(deltaz_l)
y_u = CST(data_upper[:, 0], 1, deltasz=deltaz_u, Au=Au)
y_l = CST(data_lower[:, 0], 1, deltasz=deltaz_l, Al=Al)
plt.figure()
plt.scatter(data_upper[:, 0], data_upper[:, 1], label='raw_upper')
plt.scatter(data_lower[:, 0], data_lower[:, 1], label='raw_lower')
plt.plot(data_upper[:, 0], y_u, label='upper')
plt.plot(data_lower[:, 0], y_l, label='lower')
plt.legend()
plt.show()
pickle.dump(f, open("single_airfoil_fit.p", "wb"))

for i in range(5):
    print('Au%i' % i, min(np.array(f['Au'])[:, i]),
          max(np.array(f['Au'])[:, i]))
for i in range(5):
    print('Al%i' % i, min(np.array(f['Al'])[:, i]),
          max(np.array(f['Al'])[:, i]))
