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


airfoil_database = pickle.load(open('./fitting.p', 'rb'))

# list of strings
Al_database = np.array(airfoil_database['Al'])
Au_database = np.array(airfoil_database['Au'])
dl_database = np.array(airfoil_database['dl'])
du_database = np.array(airfoil_database['du'])

airfoil = 'from_database'
altitude = 10000
chord = 1
n = 10
velocities = np.linspace(20, 65, n)
AOAs = np.linspace(0, 12, n)
AOAs, velocities = np.meshgrid(AOAs, velocities)

data = {'Names':airfoil_database['names'], 'AOA':AOAs, 'V':velocities, 'L/D':[]}

for j in range(len(Au_database)):
    print(j, airfoil_database['names'][j])
    Au = Au_database[j, :]
    Al = Al_database[j, :]
    x = create_x(1., distribution = 'linear')
    y = CST(x, chord, deltasz=[du_database[j], dl_database[j]],
                     Al=Al, Au=Au)

    xf.create_input(x, y['u'], y['l'], airfoil, different_x_upper_lower = False)
    for i in range(n):
        for k in range(n):
            AOA = AOAs[i][k]
            V = velocities[i][k]
            try:
                Data = xf.find_coefficients(airfoil, AOA,
                                            Reynolds=Reynolds(10000, V, chord),
                                            iteration=100, NACA=False,
                                            delete=True)
                lift_drag_ratio = Data['CL']/Data['CD']
            except:
                lift_drag_ratio = None
                increment = 0.01
                conv_counter = 0
                while lift_drag_ratio is None and conv_counter <2:
                    try:
                        print(increment)
                        Data_f = xf.find_coefficients(airfoil, AOA+increment,
                                                      Reynolds=Reynolds(10000, V, chord),
                                                      iteration=100, NACA=False,
                                                      delete=True)
                        Data_b = xf.find_coefficients(airfoil, AOA-increment,
                                                      Reynolds=Reynolds(10000, V, chord),
                                                      iteration=100, NACA=False,
                                                      delete=True)
                        lift_drag_ratio = .5*(Data_f['CL']/Data_f['CD'] +
                                              Data_b['CL']/Data_b['CD'])
                    except:
                        increment += increment
                    conv_counter += 1
            print(airfoil_database['names'][j], AOA, V, lift_drag_ratio)
            data['L/D'].append(lift_drag_ratio)
f = open('aerodynamics.p', 'wb')
pickle.dump(data,f)
f.close()
