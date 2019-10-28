import pickle
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

from weather.scraper.flight_conditions import properties, Airframe

# Define object
fuel = 56*6.01*0.4535
initial_mass = 1111
final_mass = initial_mass-fuel
C172_props = properties({'Cl_alpha': 5.143, 'Cl_0': 0.31,
                         'planform': 16.1651, 'density': 0.770488088,
                         'mass_min': final_mass, 'mass_max': initial_mass,
                         'incidence': 0.})
C172 = Airframe(airframe='C172', timestamp=1549036800,
                filepath='../../../weather/data/flight_plan/v_aoa_pickles/icao24s_',
                properties=C172_props)
C172.retrieve_data(load_data=True)
C172.train_pdf(1000)

# Calculating total probability
xgrid = np.linspace(-5, 35, 1000)
ygrid = np.linspace(20, 75, 1000)
X, Y = np.meshgrid(xgrid, ygrid)
Z = np.exp(C172.pdf.score_samples(np.array([X.ravel(), Y.ravel()]).T))
Z = np.reshape(Z, X.shape)
total_list = []
for i in range(len(Z)):
    # print('X', X[i])
    # print('Y', Y[:, 0])
    # print('Z', Z[i, :])
    numerator = simps(Z[i, :], X[i])
    total_list.append(numerator)
total = simps(total_list, Y[:, 0])
print('Probability total', total)

# Plot histograms
parameters = []
for i in range(200):
    sample = C172.pdf.sample(1)
    while sample[0][0]<0 or sample[0][0]>12:
        sample = C172.pdf.sample(1)
    parameters.append(sample)
C172.samples = np.array(parameters)

f = open('c172.p', 'wb')
pickle.dump(C172, f)
f.close()

# Plot PDF
C172.plot_pdf()
x, y = C172.samples.T
plt.scatter(x, y, c='k')
plt.show()
