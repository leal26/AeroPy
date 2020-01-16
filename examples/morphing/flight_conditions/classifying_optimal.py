import numpy as np
import matplotlib.pyplot as plt
from aeropy.geometry.airfoil import CST
from aeropy.morphing.camber_2D import *

import pandas
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min

# testing = 'structurally_consistent'

inverted = False
morphing_direction = 'forwards'


# ==============================================================================
# Inputs
# ==============================================================================
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

data = pandas.read_csv('optimal_map.csv')

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
# points = []
# designs = []
# for i in range(len(data.values[:, 0])):
#     print(i)
#     AC_u = list(data.values[i, 0:4])
#     Au_C, Al_C, c_C, spar_thicknesses = calculate_dependent_shape_coefficients(
#                                 AC_u,
#                                 psi_spars, Au_P, Al_P,
#                                 deltaz, c_P, morphing=morphing_direction)
#     x = np.linspace(0, c_C, 100)
#     y = CST(x, c_C, deltasz=[deltaz/2., deltaz/2.], Al=Al_C, Au=Au_C)
#     points.append(list(x) + list(y['l']) + list(y['u']))
#     designs.append({'x': x, 'yl': y['l'], 'yu': y['u']})
# points = np.array(points)
# f = open('design_optimal.p', 'wb')
# pickle.dump(designs, f)
# f.close()
# f = open('points_optimal.p', 'wb')
# pickle.dump(points, f)
# f.close()
n_clusters = 4
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
# colors = ['0.3', '0.5', '0.7']
colors = ['b', 'g', 'r', 'm', 'c']
plt.figure()
for ii in range(n_clusters):
    i = y_km[closest[ii]]
    plt.scatter(data[y_km == i, -3], data[y_km == i, -2], s=100, c=colors[ii], label=ii)
plt.scatter(data[closest, -3], data[closest, -2], s=100, c='k',
            marker='s', label='Centers')
plt.xlabel(r'Angle of Attack (${}^\circ$)')
plt.ylabel('Velocity (m/s)')
plt.legend()
plt.show()

# ==============================================================================
#  Plot results
# ==============================================================================
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


for ii in range(len(closest)):
    i = closest[ii]
    d = designs[i]
    plt.figure()

    AC_u = list(data[i, 0:4])
    print(i, data[i])
    Au_C, Al_C, c_C, spar_thicknesses = calculate_dependent_shape_coefficients(
        AC_u,
        psi_spars, Au_P, Al_P,
        deltaz, c_P, morphing=morphing_direction)
    x_c = np.linspace(0, c_C, 1000)
    y_c = CST(x_c, c_C, deltasz=[deltaz/2., deltaz/2.], Al=Al_C, Au=Au_C)

    plt.plot(d['x'], d['yu'], colors[ii], label='Children', lw=2)
    plt.plot(d['x'], d['yl'], colors[ii], label=None, lw=2)
    plt.plot(x_p, y_p['u'], 'k--', label='Parent', lw=2)
    plt.plot(x_p, y_p['l'], 'k--', label=None, lw=2)

    if morphing_direction == 'forwards':
        psi_flats = []
        intersections_x_children = [0]
        intersections_y_children = [0]
        intersections_x_parent = [0]
        intersections_y_parent = [0]
        for j in range(len(psi_spars)):
            psi_parent_j = psi_spars[j]
            # Calculate psi at landing
            # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
            psi_children_j = calculate_psi_goal(psi_parent_j, Au_P, Au_C, deltaz, c_P, c_C)
            x_children_j = psi_children_j*c_C

            # Calculate xi at landing
            temp = CST(x_children_j, c_C, [deltaz/2., deltaz/2.], Al=Al_C, Au=Au_C)
            y_children_j = temp['u']

            s = calculate_spar_direction(psi_spars[j], Au_P, Au_C, deltaz, c_C)

            # Print spars for children
            if not inverted:
                plt.plot([x_children_j, x_children_j - spar_thicknesses[j]*s[0]], [y_children_j,
                                                                                   y_children_j - spar_thicknesses[j]*s[1]], c=colors[ii], lw=2, label=None)
            else:
                plt.plot([x_children_j, x_children_j - spar_thicknesses[j]*s[0]],
                         [-y_children_j, -y_children_j + spar_thicknesses[j]*s[1]], c=colors[ii], lw=2, label=None)
            psi_flats.append(x_children_j - spar_thicknesses[j]*s[0])
            y = CST(np.array([psi_parent_j*c_P]), c_P,
                    deltasz=[deltaz/2., deltaz/2.], Al=Al_P, Au=Au_P)

            intersections_x_children.append(x_children_j - spar_thicknesses[j]*s[0])
            intersections_y_children.append(y_children_j - spar_thicknesses[j]*s[1])

            # Print spars for parents
            if not inverted:
                plt.plot([psi_parent_j*c_P, psi_parent_j*c_P],
                         [y['u'], y['u']-spar_thicknesses[j]], 'k--', lw=2, label=None)
            else:
                plt.plot([psi_parent_j*c_P, psi_parent_j*c_P], [-y['u'], -
                                                                y['u']+spar_thicknesses[j]], 'k--', lw=2, label=None)

            intersections_x_parent.append(psi_parent_j*c_P)
            intersections_y_parent.append(y['u']-spar_thicknesses[j])
    elif morphing_direction == 'backwards':
        # For backwards, goal is the parent and deformed is children
        for i in range(len(psi_spars)):
            psi_i = psi_spars[i]
            # Calculate psi at landing
            psi_goal_i = calculate_psi_goal(psi_i, Au_C, Au_P, deltaz, c_C, c_P)
            x_goal_i = psi_goal_i*c_P
            # Calculate xi at landing
            temp = CST(x_goal_i, c_P, [deltaz/2., deltaz/2.], Al=Al_P, Au=Au_P)
            y_goal_i = temp['u']

            # calculate spar direction
            s = calculate_spar_direction(psi_i, Au_C, Au_P, deltaz, c_P)

            plt.plot([x_goal_i, x_goal_i - spar_thicknesses[i]*s[0]],
                     [y_goal_i, y_goal_i - spar_thicknesses[i]*s[1]], 'k--')

            y = CST(np.array([psi_i*c_C]), c_C, deltasz=[deltaz/2., deltaz/2.], Al=Al_C, Au=Au_C)

            plt.plot([psi_i*c_C, psi_i*c_C], [y['u'], y['u'] -
                                              spar_thicknesses[i]], colors[ii], lw=2, label=None)

    plt.xlabel('$\psi^p$', fontsize=14)
    plt.ylabel(r'$\zeta^p$', fontsize=14)
    plt.ylim([-0.06, 0.17])
    plt.grid()
    plt.axis('off')
    plt.gca().set_aspect('equal', adjustable='box')
    # plt.legend(loc=1)
    plt.show()

# if morphing_direction == 'forwards':
#     print('chords', c_P, c_C)
#     # Calculate initial lengths
#     strains, av_strains = calculate_strains(Au_P, Al_P, c_P, Au_C, Al_C, c_C, deltaz, psi_spars)
#
#     intersections_x_children.append(c_C)
#     intersections_y_children.append(0)
#     intersections_x_parent.append(c_P)
#     intersections_y_parent.append(0)
#     # Wire lengths
#     for i in range(len(intersections_x_children)-1):
#         length_parent = math.sqrt((intersections_x_parent[i]-intersections_x_parent[i+1])**2 +
#                                   (intersections_y_parent[i]-intersections_y_parent[i+1])**2)
#         length_children = math.sqrt((intersections_x_children[i]-intersections_x_children[i+1])**2 +
#                                     (intersections_y_children[i]-intersections_y_children[i+1])**2)
#         print((length_children-length_parent)/length_parent)
