import aeropy.xfoil_module as xf
from aeropy.geometry.airfoil import CST, create_x
from aeropy.morphing.camber_2D import *
from aeropy.aero_module import air_properties, Reynolds, LLT_calculator

import numpy as np
import matplotlib.pyplot as plt
import pandas
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
import scipy


def coefficient_LLT(AC, velocity, AOA):
    Au_P = [0.1828, 0.1179, 0.2079, 0.0850, 0.1874]
    Al_P = Au_P
    deltaz = 0

    # Determine children shape coeffcients
    AC_u = list(data.values[i, 0:4])
    Au_C, Al_C, c_C, spar_thicknesses = calculate_dependent_shape_coefficients(
        AC_u,
        psi_spars, Au_P, Al_P,
        deltaz, c_P, morphing=morphing_direction)

    # Calculate aerodynamics for that airfoil
    airfoil = 'optimal'
    x = create_x(1., distribution='linear')
    y = CST(x, 1., [deltaz/2., deltaz/2.], Al=Al_C, Au=Au_C)
    # Create file for Xfoil to read coordinates
    xf.create_input(x, y['u'], y['l'], airfoil, different_x_upper_lower=False)
    Data = xf.find_coefficients(airfoil, AOA, Reynolds=Reynolds(
        10000, velocity, c_C), iteration=100, NACA=False)
    deviation = 0.001
    while Data['CL'] is None:
        Data_aft = xf.find_coefficients(
            airfoil, AOA*deviation, Reynolds=Reynolds(10000, velocity, c_C), iteration=100, NACA=False)
        Data_fwd = xf.find_coefficients(
            airfoil, AOA*(1-deviation), Reynolds=Reynolds(10000, velocity, c_C), iteration=100, NACA=False)
        try:
            for key in Data:
                Data[key] = (Data_aft[key] + Data_fwd[key])/2.
        except:
            deviation += deviation
    alpha_L_0 = xf.find_alpha_L_0(airfoil, Reynolds=0, iteration=100, NACA=False)

    coefficients = LLT_calculator(alpha_L_0, Data['CD'], N=100, b=span, taper=1.,
                                  chord_root=chord_root, alpha_root=AOA, V=velocity)
    lift = coefficients['C_L']
    drag = coefficients['C_D']

    return lift, drag


def aircraft_range_varying_AOA(data, lift, velocity):
    aoa, v, lift_to_drag = data
    f_LD = scipy.interpolate.interp2d(aoa, v, lift_to_drag, kind='cubic')
    f_L = scipy.interpolate.interp2d(aoa, v, lift, kind='cubic')
    def to_integrate(weight):
        # velocity = 0.514444*108 # m/s (113 KTAS)
        def calculate_AOA(velocity):
            def residual(AOA):
                CL = f_L(AOA, velocity)
                span = 11
                AR = 7.32
                chord_root = span/AR
                dynamic_pressure = .5*density*velocity**2*(span*chord_root)
                return abs(CL - weight/dynamic_pressure)
            res = scipy.optimize.minimize(residual, 0)
            return res.x

        AOA = calculate_AOA(velocity)
        lift_to_drag = f_LD(AOA, velocity)

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

    g = 9.81  # kg/ms
    fuel = 56*6.01*0.4535*g
    initial_weight = 1111*g
    final_weight = initial_weight-fuel
    return scipy.integrate.quad(to_integrate, final_weight, initial_weight)[0]


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

f = open('wing.p', 'rb')
data = pickle.load(f)
f.close()

# print(data)
# def f_LD(AOA, velocity):
# return scipy.interpolate.griddata(data.values[:,-5:-3], data.values[i,-1], np.array([[AOA, velocity],]), method='cubic')
# def f_L(AOA, velocity):
# print(AOA, velocity)
# print(data.values[:,-5:-3])
# print(np.array([[AOA, velocity],]))
# return scipy.interpolate.griddata(data.values[:,-5:-3], data.values[i,-3], np.array([[AOA, velocity],]), method='cubic')

C172 = pickle.load(open('c172_2.p', 'rb'))

data = {}
states = ['nonmorphed', 'morphed']
concepts = ['UAG8814320', 'NACA0012', 'NACA4415', 'NACA641212']
for state in ['nonmorphed', 'morphed']:
    data[state] = {}
    for concept in concepts:
        if state == 'nonmorphed':
            mat =  scipy.io.loadmat('./'+state+'/'+concept+'.mat')
            print(mat.keys())
            # AOA, V = np.meshgrid(mat['aoa'], mat['V'])
            # LD_ratio = mat['lift_to_drag']
            # data[state][concept] = [AOA, V, LD_ratio]
            # print(state, concept, max(LD_ratio.ravel()), expected_MonteCarlos(data[state][concept], C172)[0])
        elif state == 'morphed':
            if concept == 'UAG8814320':
                mat =  scipy.io.loadmat('./'+state+'/comparisons.mat')['glider'][0][0]
            else:
                mat =  scipy.io.loadmat('./'+state+'/comparisons.mat')[concept][0][0]
            AOA = mat[-8].T
            V = mat[-7].T
            LD_ratio = mat[-2].T
            LIFT = mat[-1].T
            print(LIFT)

            data[state][concept] = [AOA, V, LD_ratio]

            ranges = []
            velocity = np.linspace(20, 65, 100)
            for i in range(len(velocity)):
                range_i = aircraft_range_varying_AOA(data[state][concept], LIFT, velocity[i])
                print(state, concept, i, velocity[i], range_i)
                ranges.append(range_i)

plt.figure()
plt.plot(velocity, ranges, 'k', lw=2)
plt.xlim(min(velocity), max(velocity))
plt.ylim(min(ranges), max(ranges))
plt.xlabel('Velocity (m/s)')
plt.ylabel('Range (km)')
plt.show()
