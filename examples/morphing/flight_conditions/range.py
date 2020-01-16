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

def aircraft_range(cl, cd, velocity):
    # velocity = 0.514444*108 # m/s (113 KTAS)
    lift_to_drag = cl/cd
    g = 9.81 # kg/ms
    span = 10.9728
    RPM = 1800
    a = 0.3089 # (lb/hr)/BTU
    b = 0.008*RPM+19.607 # lb/hr 

    lbhr_to_kgs = 0.000125998
    BHP_to_watt = 745.7

    eta = 0.85

    initial_weight = 1111
    fuel = 56*6.01*0.4535
    final_weight = initial_weight-fuel
    elliptical_area = 2*np.sqrt(span/2)*3.1415/4.
    drag = .5*velocity**2*elliptical_area*cd
    thrust=drag

    power_SI = thrust*velocity/eta
    power_BHP = power_SI/BHP_to_watt
    mass_flow = (a*power_BHP + b)
    mass_flow_SI = mass_flow*lbhr_to_kgs

    SFC = mass_flow_SI/thrust
    R = velocity/g/SFC*lift_to_drag*np.log(initial_weight/final_weight)
    return R*0.0005399
    
def aircraft_range_NC(cl, cd, velocity):
    def to_integrate(weight):
        # velocity = 0.514444*108 # m/s (113 KTAS)
        lift_to_drag = cl/cd

        span = 10.9728
        RPM = 1800
        a = 0.3089 # (lb/hr)/BTU
        b = 0.008*RPM+19.607 # lb/hr 

        lbhr_to_kgs = 0.000125998
        BHP_to_watt = 745.7

        eta = 0.85

        thrust=weight/lift_to_drag

        power_SI = thrust*velocity/eta
        power_BHP = power_SI/BHP_to_watt
        mass_flow = (a*power_BHP + b)
        mass_flow_SI = mass_flow*lbhr_to_kgs

        SFC = mass_flow_SI/thrust
        dR = velocity/g/SFC*lift_to_drag/weight
        print(power_BHP)
        return dR*0.0005399
    g = 9.81 # kg/ms
    fuel = 56*6.01*0.4535*g
    initial_weight = 1111*g
    final_weight = initial_weight-fuel
    return scipy.integrate.quad(to_integrate, final_weight, initial_weight)[0]
inverted = False
morphing_direction = 'forwards'

def aircraft_range_LLT(AC, velocity, AOA):
    def to_integrate(weight):
        # velocity = 0.514444*108 # m/s (113 KTAS)

        span = 10.9728
        RPM = 1800
        a = 0.3089 # (lb/hr)/BTU
        b = 0.008*RPM+19.607 # lb/hr 

        lbhr_to_kgs = 0.000125998
        BHP_to_watt = 745.7

        eta = 0.85

        thrust=weight/lift_to_drag

        power_SI = thrust*velocity/eta
        power_BHP = power_SI/BHP_to_watt
        mass_flow = (a*power_BHP + b)
        mass_flow_SI = mass_flow*lbhr_to_kgs

        SFC = mass_flow_SI/thrust
        dR = velocity/g/SFC*lift_to_drag/weight
        return dR*0.001 #*0.0005399
        

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
    x = create_x(1., distribution = 'linear')
    y = CST(x, 1., [deltaz/2., deltaz/2.], Al=Al_C, Au=Au_C)
    # Create file for Xfoil to read coordinates
    xf.create_input(x, y['u'], y['l'], airfoil, different_x_upper_lower = False)
    Data = xf.find_coefficients(airfoil, AOA, Reynolds=Reynolds(10000, velocity, c_C), iteration=100, NACA=False)
    deviation = 0.001
    while Data['CL'] is None:
        Data_aft = xf.find_coefficients(airfoil, AOA*deviation, Reynolds=Reynolds(10000, velocity, c_C), iteration=100, NACA=False)
        Data_fwd = xf.find_coefficients(airfoil, AOA*(1-deviation), Reynolds=Reynolds(10000, velocity, c_C), iteration=100, NACA=False)
        try:
            for key in Data:
                Data[key] = (Data_aft[key] + Data_fwd[key])/2.
        except:
            deviation += deviation
    alpha_L_0 = xf.find_alpha_L_0(airfoil, Reynolds=0, iteration=100, NACA=False)
    
    coefficients = LLT_calculator(alpha_L_0, Data['CD'], N=100, b=span, taper=1.,
                   chord_root=chord_root, alpha_root=AOA, V=velocity)
    lift_to_drag = coefficients['C_L']/coefficients['C_D']
    
    g = 9.81 # kg/ms
    fuel = 56*6.01*0.4535*g
    initial_weight = 1111*g
    final_weight = initial_weight-fuel
    return scipy.integrate.quad(to_integrate, final_weight, initial_weight)[0]
    
# ==============================================================================
# Inputs
# ==============================================================================
altitude = 10000 # ft
air_props = air_properties(altitude, unit='feet')
density = air_props['Density']

data = pandas.read_csv('performance_grid.csv')
psi_spars = [0.1, 0.3, 0.6, 0.8]
c_P = 1.0
span = 11
chord_root = span/16.2
ranges = []
for i in range(len(data.values)):
    AC = data.values[i,0:4]
    velocity = data.values[i,-4]
    AOA = data.values[i,-5]
    range_i = aircraft_range_LLT(AC, velocity, AOA)
    ranges.append(range_i)
    print(i, range_i)
data['Range'] = ranges

x = data['AOA']
y = data['V']
z = data['Range']
N = 100
xi = np.linspace(x.min(), x.max(), N)
yi = np.linspace(y.min(), y.max(), N)
zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

fig = plt.figure()
plt.contourf(xi, yi, zi)
plt.xlabel("Angle of attack ($^{\circ}$)")
plt.ylabel("Velocity (m/s)")
plt.colorbar(label='Range (km)')
plt.show()
