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
    try:
        pdf = airFrame.pdf.score_samples(np.vstack([alpha.ravel(), velocity.ravel()]).T)
    except:
        print(alpha)
        print(velocity)
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

C172 = pickle.load(open('C172_new.p', 'rb'))
data = pickle.load(open('aerodynamics_3.p', 'rb'))
airfoil_database = pickle.load(open('fitting.p', 'rb'))

# list of strings
Al_database = np.array(airfoil_database['Al'])
Au_database = np.array(airfoil_database['Au'])
dl_database = np.array(airfoil_database['dl'])
du_database = np.array(airfoil_database['du'])

airfoil = 'from_database_3'
altitude = 10000
chord = 1

[AOAs, velocities] = C172.samples.T
AOAs = AOAs[0]
velocities = velocities[0]

# data = {'Names':airfoil_database['names'], 'AOA':AOAs, 'V':velocities,
#         'L/D':[], 'Expected':[]}
print(len(data['L/D']))

for j in range(len(data['L/D'])):
    if len(data['L/D'][j]) == len(AOAs.flatten()):
        data_i = np.array([AOAs.flatten(), velocities.flatten(),
                           np.array(data['L/D'][j])])
        data['Expected'].append(expected(data_i, C172))
    else:
        data['Expected'].append(0)

max_value = max(data['Expected'])
max_index = data['Expected'].index(max_value)
# df = pd.DataFrame(np.array([data['Names'], data['Expected']]).T,
#                   columns = ['Names',  'Expected'])
df = pd.DataFrame(np.array([data['Expected']]).T,
                  columns = ['Expected'])
df['Expected'] = df['Expected'].astype(float)
df.hist(column = 'Expected', cumulative=True, normed=True, bins=20)
plt.show()
print(max_index )
print(df.max())
print(airfoil_database['names'][max_index-3:max_index+3])
# print(data['L/D'][max_index])

C172.plot_pdf()
x, y = C172.samples.T
plt.scatter(x, y, c='k')
plt.show()
