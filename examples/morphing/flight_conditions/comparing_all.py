
from aeropy.geometry.airfoil import CST
from aeropy.morphing.camber_2D import *

import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import interpolate
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min

# # testing = 'structurally_consistent'

# inverted = False
# morphing_direction = 'forwards'


# airfoil_database = pickle.load(open('../2D/fitting.p', 'rb'))

# # list of strings
# Al_database = np.array(airfoil_database['Al'])
# Au_database = np.array(airfoil_database['Au'])
# dl_database = np.array(airfoil_database['dl'])
# du_database = np.array(airfoil_database['du'])

# # Parameter
# c_P = 1.  # m
# deltaz = 0.*c_P  # m

# # Avian wing, order 5
# # Au_P = [0.23993240191629417, 0.34468227138908186, 0.18125405377549103,
# #         0.35371349126072665, 0.2440815012119143, 0.25724974995738387]
# # Al_P = [0.18889012559339036, -0.24686758992053115, 0.077569769493868401,
# #         -0.547827192265256, -0.0047342206759065641, -0.23994805474814629]
# # NACA0012
# Au_P = [0.1828, 0.1179, 0.2079, 0.0850, 0.1874]
# Al_P = Au_P
# # Passive shape coefficients for parent
# # Au_P = [.5,.4,.3]
# # Active shape coefficients for parent
# # Al_P = [.5,.1,.1]

# n = len(Au_P) - 1

# if inverted:
    # temp = Au_P
    # Au_P = list(-np.array(Al_P))
    # Al_P = list(-np.array(temp))

# # Passive shape coefficients for child

# data = pd.read_csv('optimal_map.csv')
# # Spar position for cruise (adiminesional because the chord will still be calculated)
# psi_spars = [0.1, 0.3, 0.6, 0.8]

# # ==============================================================================
# # Calculate dependent coefficients
# # ==============================================================================
# import pickle
# f = open('design_optimal.p', 'rb')
# designs = pickle.load(f)
# f.close()
# f = open('points_optimal.p', 'rb')
# points = pickle.load(f)
# f.close()

# # Find points from database that are closest to centers
# closest_error = []
# closest_name = []

# for i in range(len(data.values)):
    # x = designs[i]['x']
    # yl_morphing = designs[i]['yl']
    # yu_morphing = designs[i]['yu']
    # camber_morphing = (yu_morphing + yl_morphing)/2.
    # chord = max(x)
    # current_rmse = 1e10
    # for j in range(len(Au_database)):
        # Au = Au_database[j, :]
        # Al = Al_database[j, :]
        # y_database = CST(x, chord, deltasz=[du_database[j], dl_database[j]],
                         # Al=Al, Au=Au)

        # rmse = np.sqrt(np.sum((yl_morphing-y_database['l'])**2 +
                              # (yu_morphing-y_database['u'])**2)/(2*len(x)))

        # if rmse <= current_rmse:
            # closest_name_i = airfoil_database['names'][j]
            # closest_error_i = rmse
            # current_rmse = rmse
    # print(i, closest_name_i, closest_error_i)
    # closest_error.append(closest_error_i)
    # closest_name.append(closest_name_i)

# data['Distance'] = closest_error
# data['Name'] = closest_name

# data.to_pickle('./all_comparison.p')

def name_to_value(name, names):
    xs = len(names)*np.linspace(0, 1, len(names)+1)
    if name in names:
        return xs[list(names).index(name)]
    else:
        return xs[-1]

limit = 5

data = pd.read_pickle("./all_comparison.p")
vc = data.Name.value_counts(normalize=True)*100
vc[vc>=limit].plot(kind='bar', figsize=(6, 3), rot=0)
plt.xlabel('Closet airfoils')
plt.ylabel('Airfoil probability (%)')
x = data['AOA'].values
y = data['V'].values
z = data['Distance'].values
data_names = data['Name'].values
points = np.array([x,y]).T


X = np.linspace(min(x), max(x), 100)
Y = np.linspace(min(y), max(y), 100)
X, Y = np.meshgrid(X, Y)
points_plot = np.array([X.flatten(), Y.flatten()]).T

plt.figure()
Z = interpolate.griddata(points, z, points_plot, method='cubic')
Z = Z.reshape(X.shape)
cs = plt.contourf(X, Y, Z)

cs.cmap.set_under('k')
plt.colorbar(cs, extend='min', label = 'Eucledian distance from closest existing airfoil')

# tag = [0, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.0055]
cmap = plt.get_cmap('gist_rainbow')
# cmaplist = [cmap(i) for i in range(cmap.N)]
# cmaplist = [(0, 0, 0, 1),].append(cmaplist)
#
# # create the new map
# cmap = mpl.colors.LinearSegmentedColormap.from_list(
#     'Custom cmap', cmaplist, cmap.N)
#
# # define the bins and normalize
# bounds = np.linspace(0, 20, 21)
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#
# # make the scatter
# scat = ax.scatter(x, y, c=tag, s=np.random.randint(100, 500, 20),
#                   cmap=cmap, norm=norm)

plt.figure()
vc_above = vc[vc>=limit]
names_above = vc_above.index.values
count = vc_above.values
colors = cmap(np.linspace(0, 1, len(names_above)+1))
for i, (name, color) in enumerate(zip(names_above, colors), 1):
    print(i, name, color)
    plt.scatter(x[data_names == name], y[data_names == name], label=name, c=np.array(color))

vc_below = vc[vc<limit]
names = vc_below.index.values
for j, name in enumerate(names, 1):
    if j == 1:
        plt.scatter(x[data_names == name], y[data_names == name], c=colors[-1], label='Misc')
    else:
        plt.scatter(x[data_names == name], y[data_names == name], c=colors[-1])
plt.legend()

plt.figure()
X = np.linspace(min(x), max(x), 26)
Y = np.linspace(min(y), max(y), 30)
X, Y = np.meshgrid(X, Y)
points_plot = np.array([X.flatten(), Y.flatten()]).T
z = []
for i in range(len(x)):
    z.append(name_to_value(data_names[i], names_above))
Z = interpolate.griddata(points, z, points_plot, method='cubic')
Z = Z.reshape(X.shape)
cs = plt.contourf(X, Y, Z)
cs.cmap.set_under('k')
plt.colorbar(label = 'Airfoils')
plt.show()
