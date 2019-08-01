import numpy as np
import matplotlib.pyplot as plt
from aeropy.geometry.airfoil import CST
from aeropy.morphing.camber_2D import *

import pickle
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min

# testing = 'structurally_consistent'

inverted = False
morphing_direction = 'forwards'


airfoil_database = pickle.load(open('../2D/fitting.p', 'rb'))

# list of strings
Al_database = np.array(airfoil_database['Al'])
Au_database = np.array(airfoil_database['Au'])
dl_database = np.array(airfoil_database['dl'])
du_database = np.array(airfoil_database['du'])

# Parameter
c_P = 1.  # m
deltaz = 0.*c_P  # m

# Avian wing, order 5
# Au_P = [0.23993240191629417, 0.34468227138908186, 0.18125405377549103,
#         0.35371349126072665, 0.2440815012119143, 0.25724974995738387]
# Al_P = [0.18889012559339036, -0.24686758992053115, 0.077569769493868401,
#         -0.547827192265256, -0.0047342206759065641, -0.23994805474814629]
# NACA0012
Au_P = [0.1828, 0.1179, 0.2079, 0.0850, 0.1874]
Al_P = Au_P
# Passive shape coefficients for parent
# Au_P = [.5,.4,.3]
# Active shape coefficients for parent
# Al_P = [.5,.1,.1]

n = len(Au_P) - 1

if inverted:
    temp = Au_P
    Au_P = list(-np.array(Al_P))
    Al_P = list(-np.array(temp))

# Passive shape coefficients for child

data = pd.read_csv('optimal_map.csv')
# Spar position for cruise (adiminesional because the chord will still be calculated)
psi_spars = [0.1, 0.3, 0.6, 0.8]

# ==============================================================================
# Calculate dependent coefficients
# ==============================================================================
import pickle
f = open('design_optimal.p', 'rb')
designs = pickle.load(f)
f.close()
f = open('points_optimal.p', 'rb')
points = pickle.load(f)
f.close()

n_clusters = 3
# create kmeans object
kmeans = KMeans(n_clusters=n_clusters)

# fit kmeans object to data
kmeans.fit(points)

# save new clusters for chart
y_km = kmeans.fit_predict(points)

# Find designs closest to cluster
data = data.values
closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, points)
closest = closest[np.argsort(data[closest, -3])]

# Find points from database that are closest to centers
closest_database = []
print(closest)
for ii in range(n_clusters):
    error = 0
    i = closest[ii]
    print(i)
    x = designs[i]['x']
    yl_morphing = designs[i]['yl']
    yu_morphing = designs[i]['yu']
    camber_morphing = (yu_morphing + yl_morphing)/2.
    chord = max(x)
    current_rmse = 1e10
    jjs = np.arange(0, len(Au_database))
    # jjs = np.delete(jjs, 991)
    for j in jjs:
        Au = Au_database[j, :]
        Al = Al_database[j, :]
        y_database = CST(x, chord, deltasz=[du_database[j], dl_database[j]],
                         Al=Al, Au=Au)
        camber_database = (y_database['u'] + y_database['l'])/2.
        # rmse = np.sqrt(np.mean((camber_morphing - camber_database)**2))
        rmse = np.sqrt(np.sum((yl_morphing-y_database['l'])**2 +
                              (yu_morphing-y_database['u'])**2)/(2*len(x)))
        error += 1
        if rmse <= current_rmse:
            closest_database_i = {'x': x, 'yl': y_database['l'],
                                  'yu': y_database['u'],
                                  'name': airfoil_database['names'][j], 'index': j, 'rmse': rmse}
            current_rmse = rmse
    print(closest_database_i['name'], closest_database_i['index'],
          len(Au_database), closest_database_i['rmse'], error)
    closest_database.append(closest_database_i)

# ==============================================================================
#  Plot results
# ==============================================================================
colors = ['b', 'g', 'r', 'm', 'c']

plt.figure()
np.set_printoptions(precision=20)
x_p = np.linspace(0, c_P, 100000)
y_p = CST(x_p, c_P, deltasz=[deltaz/2., deltaz/2.], Al=Al_P, Au=Au_P)
for ii in range(len(closest)):
    i = closest[ii]
    d = designs[i]
    plt.plot(d['x'], d['yu'], colors[ii], label='%i' % ii, lw=2)
    plt.plot(d['x'], d['yl'], colors[ii], label=None, lw=2)
plt.plot(x_p, y_p['u'], 'k--', label='Parent', lw=2)
plt.plot(x_p, y_p['l'], 'k--', label=None, lw=2)
plt.xlabel('$\psi^p$', fontsize=14)
plt.ylabel(r'$\zeta^p$', fontsize=14)
plt.ylim([-0.06, 0.17])
plt.grid()
plt.gca().set_aspect('equal', adjustable='box')
plt.legend(loc=1)
plt.show()


plt.figure()
np.set_printoptions(precision=20)
x_p = np.linspace(0, c_P, 100000)
y_p = CST(x_p, c_P, deltasz=[deltaz/2., deltaz/2.], Al=Al_P, Au=Au_P)
print(closest)
for ii in range(len(closest)):
    i = closest[ii]
    print(i, ii)
    d = designs[i]
    c = closest_database[ii]
    print(c['name'])
    plt.plot(d['x'], d['yl'], colors[ii], label='%i' % ii, lw=2)
    plt.plot(d['x'], d['yu'], colors[ii], label=None, lw=2)
    plt.plot(d['x'], (d['yu']+d['yl'])/2., colors[ii], label=None, lw=2)
    plt.plot(c['x'], c['yl'], colors[ii]+'--', label=c['name'], lw=2)
    plt.plot(c['x'], c['yu'], colors[ii]+'--', label=None, lw=2)
    plt.plot(c['x'], (c['yu']+c['yl'])/2., colors[ii]+'--', label=None, lw=2)

    plt.axis('off')
    plt.xlabel('$\psi^p$', fontsize=14)
    plt.ylabel(r'$\zeta^p$', fontsize=14)
    plt.ylim([-0.06, 0.17])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend(loc=1)
    plt.show()