from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle

def poly(x, c2, c3):
    return c2*x**2 + c3*x**3

df = pd.read_excel('curved_data.xlsx', sheet_name = 'Chen (2010)')

skip_lines = 1 # -1 of what you think it should be

plt.figure()
for i in range(int(len(df.columns)/2)):
    ii = 2*i
    x = np.array(df.values[skip_lines:, ii])
    y = np.array(df.values[skip_lines:, ii+1])
    # Remove NaN
    print(x)
    x = np.array(list(x[~np.isnan(list(x))]))
    y = np.array(list(y[~np.isnan(list(y))]))

    # Fit
    popt, pcov = curve_fit(poly, x, y)
    print(popt)
    plt.plot(x, y, 'b', label = 'Raw %i' % i)
    plt.plot(x, poly(x, *popt), 'r--', label = 'Fit %i' % i)
plt.show()

def poly(x, c2, c3, c4, c5):
    return c2*x**2 + c3*x**3 + c4*x**4 + c5*x**5

plt.figure()
abaqus_x = [-6.143523e-37, 0.049518712, 0.097414128, 0.14182274, 0.18338385, 0.22280164, 0.26078638, 0.29787692, 0.33455449, 0.3694379, 0.40472263, 0.44081238, 0.47817007, 0.51734918, 0.55898702, 0.60366577, 0.65133238, 0.70004117, 0.74728715, 0.78860354, 0.82396096, 0.85445607, 0.8812418, 0.90515256, 0.92684138, 0.94658399, 0.96465063, 0.98153377, 0.99782264]
abaqus_y = [-1.0867065e-36, 0.0069836318, 0.025579287, 0.051580817, 0.08195062, 0.11507001, 0.14982793, 0.18554173, 0.22168092, 0.25601763, 0.2899411, 0.32300428, 0.35462314, 0.38394535, 0.40963244, 0.4294889, 0.44009688, 0.43760237, 0.42023277, 0.39134949, 0.35533783, 0.31509525, 0.27228072, 0.22779004, 0.18217251, 0.13567843, 0.088506527, 0.040896866, -0.0069233831]

x = np.array(abaqus_x)
y = np.array(abaqus_y)

# Fit
popt, pcov = curve_fit(poly, x, y)
print(popt)
plt.plot(x, y, 'b', label = 'Raw')
plt.plot(x, poly(x, *popt), 'r--', label = 'Fit')
plt.show()

