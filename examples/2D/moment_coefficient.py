'''Example comparing moment coefficient from xfoil relative to an in-house
code. The in-house code allows for moment calculation around regions of aifoils
(e.g. flap)'''

import matplotlib.pyplot as plt
import numpy as np

import aeropy.aero_module as ar
import aeropy.xfoil_module as xf

# For a single angle of attack
alpha = 0.

data = xf.find_pressure_coefficients('naca0012', alpha)
C_m = ar.calculate_moment_coefficient(data['x'], data['y'], data['Cp'], alpha)
data_CM = xf.find_coefficients('naca0012', alpha, delete=True)

print('calculated:', C_m)
print('objective:', data_CM['CM'])

# FOr multiple angles of attack
Cm_xfoil = []
Cm_aeropy = []

alpha_list = np.linspace(0, 10, 11)
for alpha in alpha_list:
    alpha = float(alpha)
    data = xf.find_pressure_coefficients('naca0012', alpha, delete=True)
    Cm_aeropy.append(ar.calculate_moment_coefficient(data['x'], data['y'],
                                                     data['Cp'], alpha))
    data_CM = xf.find_coefficients('naca0012', alpha, delete=True)
    Cm_xfoil.append(data_CM['CM'])
plt.plot(alpha_list, Cm_xfoil, 'b', label='XFOIL')
plt.plot(alpha_list, Cm_aeropy, 'g', label='AeroPy')
plt.legend()
plt.xlabel("Angle of attack($^{\circ}$)")
plt.ylabel("$C_m$")
plt.show()
