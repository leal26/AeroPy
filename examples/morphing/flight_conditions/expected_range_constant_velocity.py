import pickle
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
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
    f_interpolation = interp1d(velocity, metric)

    for i in range(len(airFrame.samples)):
        data_i = C172.denormalize(np.array(airFrame.samples[i,:][0])).T
        # print(data_i)
        pdf = airFrame.pdf_velocity.score_samples([[data_i[1]],])
        pdf = np.exp(pdf)
        
        try:
            metric = f_interpolation([data_i[1]])
        except(ValueError):
            print('Error')
            print(data_i)
        # print(pdf, LD_interpolated)
        expected_value += metric[0]
        pdfs.append(pdf[0])
    total_pdf = sum(pdfs)
    # print(expected_value, total_pdf)
    expected_value = expected_value/len(airFrame.samples)
    return(expected_value)
    
C172 = pickle.load(open('c172_2.p', 'rb'))

data = {}
states = ['nonmorphed', 'morphed']
concepts = ['NACA0012', 'NACA4415', 'NACA641212', 'glider']
for state in states:
    ranges = pickle.load(open('./'+state+'/'+'ranges_velocity.p', 'rb'))
    data[state] = {}
    for concept in concepts:
        if state == 'nonmorphed':
            mat =  scipy.io.loadmat('./'+state+'/'+ state + '_' + concept)
            aoa = mat['aoa'][0]
            velocity = mat['V'][0]
            print(velocity[-1])
            cl = mat['CL'][0]
            range_i = ranges[concept]
            print(state, concept, range_i[-1], expected_value(aoa, velocity, cl, range_i, C172))
        elif state == 'morphed':
            data = np.loadtxt('./'+state + '/' + state + '_' + concept + '.txt')
            aoa = data[:,0]
            velocity = data[:,1]
            print(velocity[-1])
            cl = data[:,2]
            LD_ratio = data[:,3]
            range_i = ranges[concept]
            print(state, concept, range_i[-1], expected_value(aoa, velocity, cl, range_i, C172))