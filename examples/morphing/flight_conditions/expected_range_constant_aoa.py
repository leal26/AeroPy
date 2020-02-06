import pickle
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from sklearn.neighbors.kde import KernelDensity
from scipy.interpolate import griddata, RegularGridInterpolator, interp1d

from weather.scraper.flight_conditions import properties, Airframe

def expected_value(alpha, velocity, cl, metric, airFrame):
    # data = data[data[:,0].argsort()]
    expected_value = 0
    # V = 12*45
    V = 1
    pdfs = []
    alpha = np.sort(np.unique(alpha.ravel()).T)
    velocity = np.sort(np.unique(velocity.ravel()).T)
    # metric = np.reshape(metric, [len(alpha), len(velocity)])
    f_interpolation = interp1d(alpha, metric)

    for i in range(len(airFrame.samples)):
        data_i = C172.denormalize(np.array(airFrame.samples[i,:][0])).T
        data_i = np.array([[data_i[0]],])
        # print(data_i)
        # print(data_i)
        pdf = airFrame.pdf_aoa.score_samples(data_i)
        pdf = np.exp(pdf)
        
        try:
            metric = f_interpolation(data_i[0][0])
            # print(data_i[0][0], metric)
        except(ValueError):
            print('Error')
            print(data_i)
        # print(pdf, LD_interpolated)
        expected_value += metric
        pdfs.append(pdf)
    total_pdf = sum(pdfs)
    expected_value = expected_value/len(airFrame.samples)
    return(expected_value)
    
C172 = pickle.load(open('c172_2.p', 'rb'))
data_aoa = np.reshape(C172.denormalize(C172.database)[:,0], [len(C172.denormalize(C172.database)[:,0]), 1])
print(np.shape(data_aoa))
C172.pdf_aoa = KernelDensity(kernel='gaussian', bandwidth=0.01, algorithm='ball_tree').fit(data_aoa)

data = {}
states = ['nonmorphed', 'morphed']
concepts = ['NACA0012', 'NACA4415', 'NACA641212', 'glider']
for state in states:
    ranges = pickle.load(open('./'+state+'/'+'ranges_aoa.p', 'rb'))
    data[state] = {}
    for concept in concepts:
        if state == 'nonmorphed':
            mat =  scipy.io.loadmat('./'+state+'/'+ state + '_' + concept)
            aoa = mat['aoa'][0]
            velocity = mat['V'][0]
            cl = mat['CL'][0]
            range_i = np.array(ranges[concept])
            print(state, concept, range_i[36], max(range_i[~np.isnan(range_i)]), expected_value(aoa, velocity, cl, range_i, C172))
        elif state == 'morphed':
            data = np.loadtxt('./'+state + '/' + state + '_' + concept + '.txt')
            aoa = data[:,0]
            velocity = data[:,1]
            cl = data[:,2]
            LD_ratio = data[:,3]
            print(state, concept, range_i[36], max(range_i[~np.isnan(range_i)]), expected_value(aoa, velocity, cl, range_i, C172))
print(velocity[36])