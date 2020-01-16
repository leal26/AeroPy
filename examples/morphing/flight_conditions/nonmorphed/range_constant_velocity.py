import aeropy.xfoil_module as xf
from aeropy.geometry.airfoil import CST, create_x
from aeropy.morphing.camber_2D import *
from aeropy.aero_module import air_properties, Reynolds, LLT_calculator
from scipy.interpolate import griddata, RegularGridInterpolator

import numpy as np
import matplotlib.pyplot as plt
import pandas
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
import scipy

def aircraft_range_varying_AOA(f_L, f_LD, velocity):
    def to_integrate(weight):
        # velocity = 0.514444*108 # m/s (113 KTAS)
        def calculate_AOA(velocity):
            def residual(AOA):
                
                CL = f_L([velocity, AOA[0]])[0]
                
                span = 11
                AR = 7.32
                chord_root = span/AR
                dynamic_pressure = .5*density*velocity**2*(span*chord_root)
                return abs(CL - weight/dynamic_pressure)
            if len(AOA_list) == 0:
                x0 = 0
            else:
                x0 = AOA_list[-1]
            res = scipy.optimize.minimize(residual, x0, bounds = [[0, 12],])#, options={'ftol':1e-9})
            return res.x[0]

        AOA = calculate_AOA(velocity)
        # print(AOA)
        AOA_list.append(AOA)
        lift_to_drag = f_LD([velocity, AOA])

        span = 10.9728
        RPM = 1800
        a = 0.3089  # (lb/hr)/BTU
        b = 0.008*RPM+19.607  # lb/hr

        lbhr_to_kgs = 0.000125998
        BHP_to_watt = 745.7

        eta = 0.85
        
        thrust = weight/lift_to_drag

        power_SI = thrust*velocity/eta
        power_BHP = power_SI/BHP_to_watt
        mass_flow = (a*power_BHP + b)
        mass_flow_SI = mass_flow*lbhr_to_kgs

        SFC = mass_flow_SI/thrust
        dR = velocity/g/SFC*lift_to_drag/weight
        return dR*0.001  # *0.0005399

    AOA_list = []
    g = 9.81  # kg/ms
    fuel = 56*6.01*0.4535*g
    initial_weight = 1111*g
    final_weight = initial_weight-fuel
    x = np.linspace(final_weight, initial_weight, 100)
    y = []
    for x_i in x:
        y.append(to_integrate(x_i)[0])
    range = scipy.integrate.simps(y, x)
    return range, AOA_list


# ==============================================================================
# Inputs
# ==============================================================================
altitude = 10000  # ft
air_props = air_properties(altitude, unit='feet')
density = air_props['Density']

# data = pandas.read_csv('performance_grid.csv')
# psi_spars = [0.1, 0.3, 0.6, 0.8]
# c_P = 1.0

# ranges = []

# for i in range(len(data.values)):
# AC = data.values[i,0:4]
# velocity = data.values[i,-4]
# AOA = data.values[i,-5]
# cl= data.values[i,-3]
# cd = data.values[i,-2]
# CL, CD = coefficient_LLT(AC, velocity, AOA)
# data.values[i, -3] = CL
# data.values[i, -2] = CD
# data.values[i, -1] = CL/CD
# print(i, CL, CD)
# data = data.drop_duplicates()

import pickle
# f = open('wing.p', 'wb')
# pickle.dump(data, f)
# f.close()
state = 'nonmorphed'
concepts = ['NACA0012', 'NACA4415', 'NACA641212', 'glider']

#
# plt.figure()
# for concept in concepts:
    # mat =  scipy.io.loadmat(state + '_' + concept)
    # aoa = mat['aoa'][0]
    # velocity = mat['V'][0]
    # cl = mat['CL'].T
    # LD_ratio = mat['lift_to_drag']
    # # print(aoa)
    # # print(velocity)
    # # print(cl)
    # f_LD = RegularGridInterpolator((velocity, aoa), LD_ratio, fill_value = 0, bounds_error = False)
    # f_L = RegularGridInterpolator((velocity, aoa), cl, fill_value = 0, bounds_error = False)
    # velocity = [20]
    # aoas = np.linspace(0,12,1000)
    # for i in range(len(velocity)):
        # data_i = np.array([velocity[i]*np.ones(np.shape(aoas)), aoas]).T
        # plt.plot(aoas, f_L(data_i), label = concept)
        # # plt.scatter(aoas, f_L((aoas, velocity[i]*np.ones(np.shape(aoas)))))
    # plt.legend()
# plt.ylabel('cl')
# plt.show()

# plt.figure()
# for concept in concepts:
    # mat =  scipy.io.loadmat(state + '_' + concept)
    # aoa = mat['aoa'][0]
    # velocity = mat['V'][0]
    # cl = mat['CL'].T
    # LD_ratio = mat['lift_to_drag']

    # f_LD = RegularGridInterpolator((velocity, aoa), LD_ratio, fill_value = 0, bounds_error = False)
    # f_L = RegularGridInterpolator((velocity, aoa), cl, fill_value = 0, bounds_error = False)
    # velocity = [20]
    # aoas = np.linspace(0,12,100)
    # for i in range(len(velocity)):
        # data_i = np.array([velocity[i]*np.ones(np.shape(aoas)), aoas]).T
        # plt.plot(aoas, f_LD(data_i), label = concept)
        # # plt.scatter(aoas, f_LD((aoas, velocity[i]*np.ones(np.shape(aoas)))))
    # plt.legend()
# plt.ylabel('Lift-to-drag ratio')
# plt.show()

range_data = {}
plt.figure()
for concept in concepts:
    mat =  scipy.io.loadmat(state + '_' + concept)
    aoa = mat['aoa'][0]
    velocity = mat['V'][0]
    cl = mat['CL'].T
    LD_ratio = mat['lift_to_drag']

    f_LD = RegularGridInterpolator((velocity, aoa), LD_ratio, fill_value = 0, bounds_error = False)
    f_L = RegularGridInterpolator((velocity, aoa), cl, fill_value = 0, bounds_error = False)

    # velocity = np.linspace(20, 65, 7)
    # plt.figure()
    # aoas = np.linspace(0,12,1000)
    # for i in range(len(velocity)):
        # data_i = np.array([velocity[i]*np.ones(np.shape(aoas)), aoas]).T
        # plt.plot(aoas, f_L(data_i), label = velocity[i])
        # # plt.scatter(aoas, f_L((aoas, velocity[i]*np.ones(np.shape(aoas)))))
    # plt.legend()
    # plt.show()

    # plt.figure()
    # aoas = np.linspace(0,12,1000)
    # for i in range(len(velocity)):
        # data_i = np.array([velocity[i]*np.ones(np.shape(aoas)), aoas]).T
        # plt.plot(aoas, f_LD(data_i), label = velocity[i])
        # # plt.scatter(aoas, f_LD((aoas, velocity[i]*np.ones(np.shape(aoas)))))
    # plt.legend()
    # plt.show()

    ranges = []
    # velocity = np.linspace(20, 60, 5)
    for i in range(len(velocity)):
        range_i, AOA_i = aircraft_range_varying_AOA(f_L, f_LD, velocity[i])
        # plt.plot(np.arange(len(AOA_i)), AOA_i, label=velocity[i])
        # plt.scatter(np.arange(len(AOA_i)),AOA_i)
        print(i, velocity[i], range_i)
        ranges.append(range_i)
    # print(velocity[36])
    range_data[concept] = ranges
    plt.plot(velocity, ranges, lw=2, label=concept)
f = open('ranges_velocity.p', 'wb')
pickle.dump(range_data, f)
f.close()



# plt.xlim(min(velocity), max(velocity))
# plt.ylim(min(ranges), max(ranges))
plt.xlabel('Velocity (m/s)')
plt.ylabel('Range (km)')
plt.legend()
plt.show()
