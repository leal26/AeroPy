import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import interpolate
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min

import aeropy.xfoil_module as xf
from aeropy.aero_module import Reynolds
from aeropy.geometry.airfoil import CST, create_x

import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

from weather.scraper.flight_conditions import properties, Airframe

def expected(data, airFrame):
    # data = data[data[:,0].argsort()]
    alpha, velocity, lift_to_drag = data
    pdf = airFrame.pdf.score_samples(np.vstack([alpha.ravel(), velocity.ravel()]).T)
    pdf = np.exp(pdf.reshape(lift_to_drag.shape))
    expected_value = 0
    numerator_list = []
    denominator_list = []
    N = len(alpha.ravel())
    # V = 12*45
    V = 1
    total_pdf = sum(pdf)
    for i in range(len(lift_to_drag)):
        if lift_to_drag[i] is not None:
            expected_value += (V/total_pdf)*pdf[i]*lift_to_drag[i]
    return(expected_value)

# Define object
C172 = pickle.load(open('C172.p', 'rb'))
airfoil_database = pickle.load(open('../2D/fitting.p', 'rb'))
Al_database = np.array(airfoil_database['Al'])
Au_database = np.array(airfoil_database['Au'])
dl_database = np.array(airfoil_database['dl'])
du_database = np.array(airfoil_database['du'])

airfoil = 'from_database_3'
altitude = 10000
chord = 1


grids = [160, 180, 200]
data = {'Mean':[], 'STD':[], 'grids': grids}
j = 1105
print(j, airfoil_database['names'][j])
for n in grids:
    expected_data = []
    for k in range(7):
        # Plot histograms
        parameters = []
        for i in range(n):
            sample = C172.pdf.sample(1)
            while sample[0][0]<0 or sample[0][0]>12:
                sample = C172.pdf.sample(1)
            parameters.append(sample)
        C172.samples = np.array(parameters)
        [AOAs, velocities] = C172.samples.T
        AOAs = AOAs[0]
        velocities = velocities[0]

        # data[-1]['V'].append(velocities)
        # data[-1]['AOA'].append(AOAs)

        Au = Au_database[j, :]
        Al = Al_database[j, :]
        x = create_x(1., distribution = 'linear')
        y = CST(x, chord, deltasz=[du_database[j], dl_database[j]],
                         Al=Al, Au=Au)

        xf.create_input(x, y['u'], y['l'], airfoil, different_x_upper_lower = False)
        LDs = []
        for i in range(len(AOAs)):
            AOA = AOAs[i]
            V = velocities[i]
            try:
                Data = xf.find_coefficients(airfoil, AOA,
                                            Reynolds=Reynolds(10000, V, chord),
                                            iteration=100, NACA=False,
                                            delete=True)
                lift_drag_ratio = Data['CL']/Data['CD']
            except:
                lift_drag_ratio = None
                increment = 0.1
                conv_counter = 0
                while lift_drag_ratio is None and conv_counter <2:
                    print(increment)
                    Data_f = xf.find_coefficients(airfoil, AOA*(1+increment),
                                                  Reynolds=Reynolds(10000, V*(1+increment), chord),
                                                  iteration=100, NACA=False,
                                                  delete=True)
                    Data_b = xf.find_coefficients(airfoil, AOA*(1-increment),
                                                  Reynolds=Reynolds(10000, V*(1-increment), chord),
                                                  iteration=100, NACA=False,
                                                  delete=True)
                    print(Data_f['CL'], Data_f['CD'])
                    print(Data_b['CL'], Data_b['CD'])
                    try:
                        lift_drag_ratio = .5*(Data_f['CL']/Data_f['CD'] +
                                              Data_b['CL']/Data_b['CD'])
                    except(TypeError):
                        increment += 0.1
                    conv_counter += 1
            print(i, AOA, V, lift_drag_ratio)
            LDs.append(lift_drag_ratio)

        data_i = np.array([AOAs.flatten(), velocities.flatten(),
                           LDs])
        expected_data.append(expected(data_i, C172))
    data['Mean'].append(np.mean(np.array(expected_data)))
    data['STD'].append(np.std(np.array(expected_data)))
print(data)
df = pd.DataFrame(data)
df.to_pickle('convergence.p')
plt.figure()
plt.errorbar(grids, data['Mean'], data['STD'], linestyle='None', capsize = 2, marker = ".")
plt.xlabel('Sample size')
plt.ylabel('Expected value')
plt.show()
