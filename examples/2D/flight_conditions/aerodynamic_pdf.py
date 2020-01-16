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
    alpha, V, lift_to_drag = data

    pdf = airFrame.pdf.score_samples(np.vstack([alpha.ravel(), V.ravel()]).T)
    pdf = np.exp(pdf.reshape(lift_to_drag.shape))
    expected_value = 0
    numerator_list = []
    denominator_list = []
    for i in range(len(lift_to_drag)):
        numerator = simps(lift_to_drag[i]*pdf[i], alpha[i])
        denominator = simps(pdf[i], alpha[i])
        numerator_list.append(numerator)
        denominator_list.append(denominator)
    numerator = simps(numerator_list, V[:, 0])
    denominator = simps(denominator_list, V[:, 0])
    expected_value = numerator/denominator
    return(expected_value)

C172 = pickle.load(open('C172_new.p', 'rb'))

airfoil_database = pickle.load(open('fitting.p', 'rb'))

# list of strings
Al_database = np.array(airfoil_database['Al'])
Au_database = np.array(airfoil_database['Au'])
dl_database = np.array(airfoil_database['dl'])
du_database = np.array(airfoil_database['du'])

airfoil = 'from_database_5'
altitude = 10000
chord = 1

[AOAs, velocities] = C172.samples.T
AOAs = AOAs[0]
velocities = velocities[0]

# data = {'Names':airfoil_database['names'], 'AOA':AOAs, 'V':velocities,
#         'L/D':[], 'Expected':[]}
f = open('aerodynamics_3.p', 'rb')
data = pickle.load(f)
f.close()

for j in range(240, len(Au_database)):
    data['L/D'].append([])
    print(j, airfoil_database['names'][j])
    Au = Au_database[j, :]
    Al = Al_database[j, :]
    x = create_x(1., distribution = 'linear')
    y = CST(x, chord, deltasz=[du_database[j], dl_database[j]],
                     Al=Al, Au=Au)

    xf.create_input(x, y['u'], y['l'], airfoil, different_x_upper_lower = False)
    for i in range(len(AOAs)):
        AOA = AOAs[i]
        V = velocities[i]
        AOA, V = C172.denormalize(np.array([AOA, V]).T)
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
            while lift_drag_ratio is None and conv_counter <3:
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
        print(airfoil_database['names'][j], AOA, V, lift_drag_ratio)
        data['L/D'][-1].append(lift_drag_ratio)
        if data['L/D'][-1].count(None) > 3:
            break
    f = open('aerodynamics_4.p', 'wb')
    pickle.dump(data,f)
    f.close()
