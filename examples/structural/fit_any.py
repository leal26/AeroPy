from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle

def poly(x, c2, c3):
    return c2*x**2 + c3*x**3

df = pd.read_excel('curved_data.xlsx', sheet_name = 'FEA (3D, c3=0.25, F=-10)')

skip_lines = 0 # -1 of what you think it should be

plt.figure()
for i in range(1): #range(int(len(df.columns)/2)):
    ii = 2*i
    x = np.array(df.values[skip_lines:, ii])
    y = np.array(df.values[skip_lines:, ii+1])
    # Remove NaN
    x = np.array(list(x[~np.isnan(list(x))]))
    y = np.array(list(y[~np.isnan(list(y))]))

    # Fit
    popt, pcov = curve_fit(poly, x, y)
    print(popt)
    plt.plot(x, y, 'b', label = 'Raw %i' % i)
    plt.plot(x, poly(x, *popt), 'r--', label = 'Fit %i' % i)
plt.show()
