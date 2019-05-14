import numpy as np
import matplotlib.pyplot as plt
import math

import aeropy.xfoil_module as xf
import aeropy.AeroPy as ar
#    print find_3D_coefficients(airfoil='naca0012', alpha=1.)
alpha = 0.
x_hinge = 0.25/0.6175
deflection = -math.pi/2.  # 0.17453292519943295 #0.0010573527055

# generate original airfoil
airfoil = "naca0012"
xf.call(airfoil, output='Coordinates')
filename = xf.file_name(airfoil, output='Coordinates')
Data = xf.output_reader(filename, output='Coordinates', header=['x', 'y'])

Cm = ar.calculate_flap_moment(Data['x'], Data['y'], alpha, x_hinge, deflection)


V = 10
altitude = 10000  # feet

Reynolds = ar.Reynolds(altitude, V, 1.0)

deflection_list = list(np.linspace(5, 30, 4))
alpha_list = list(np.linspace(0, 15, 20))

# Calculate coefficients for without flap
CL_list = []
CD_list = []
ratio_list = []
CM_list = []

for alpha_i in alpha_list:
    Data_0 = xf.find_coefficients('naca0012', alpha_i, Reynolds=Reynolds,
                                  iteration=200)
    CL_list.append(Data_0['CL'])
    CD_list.append(Data_0['CD'])
    ratio_list.append(Data_0['CL']/Data_0['CD'])
    CM_list.append(Data_0['CM'])
All_data = {0: {r'$c_m$': CM_list, r'$c_l$': CL_list,
                r'$c_d$': CD_list,
                r'$c_l/c_d$': ratio_list}}  # :Data_0['CL']/Data_0['CD']}}

# Calculate foeccifient when using flap
for deflection_i in deflection_list:
    CL_list = []
    CD_list = []
    ratio_list = []
    CM_list = []
    for alpha_i in alpha_list:
        flap_data = ar.calculate_flap_coefficients(Data['x'], Data['y'],
                                                   alpha_i, x_hinge,
                                                   deflection_i,
                                                   Reynolds=Reynolds)
        CL_list.append(flap_data['CL'])
        CD_list.append(flap_data['CD'])
        ratio_list.append(flap_data['CL']/flap_data['CD'])
        CM_list.append(flap_data['CM'])
    All_data[deflection_i] = {r'$c_m$': CM_list, r'$c_l$': CL_list,
                              r'$c_d$': CD_list, r'$c_l/c_d$': ratio_list}
for key in [r'$c_m$', r'$c_l$', r'$c_d$', r'$c_l/c_d$']:
    plt.figure()
    for deflection_i in [0] + deflection_list:
        plt.plot(alpha_list, All_data[deflection_i][key],
                 label=r'$\theta$ = %.0f' % deflection_i)
    plt.legend(loc="best")
    plt.xlabel(r'$\alpha$', fontsize=22)
    plt.ylabel(key, fontsize=22)
    plt.grid()

plt.figure()
for deflection_i in [0] + deflection_list:
    plt.plot(All_data[deflection_i][r'$c_d$'],
             All_data[deflection_i][r'$c_l$'],
             label=r'$\theta$ = %.0f' % deflection_i)
plt.legend(loc="best")
plt.xlabel(r'$c_d$', fontsize=22)
plt.ylabel(r'$c_l$', fontsize=22)
plt.grid()
