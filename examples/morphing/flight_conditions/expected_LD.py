import pickle
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import griddata, RegularGridInterpolator

from weather.scraper.flight_conditions import properties, Airframe

def expected_LD(alpha, velocity, cl, lift_to_drag, airFrame):
    # data = data[data[:,0].argsort()]
    expected_value = 0
    # V = 12*45
    V = 1
    pdfs = []
    alpha = np.sort(np.unique(alpha.ravel()).T)
    velocity = np.sort(np.unique(velocity.ravel()).T)
    lift_to_drag = np.reshape(lift_to_drag, [len(alpha), len(velocity)])
    f_interpolation = RegularGridInterpolator((alpha, velocity), lift_to_drag)
    for i in range(len(airFrame.samples)):
        pdf = airFrame.pdf.score_samples(airFrame.samples[i,:])
        pdf = np.exp(pdf)
        # print(alpha.ravel().shape)
        # print(velocity.ravel().shape)
        # print(lift_to_drag.ravel().shape)
        # print(pdf)
        # print(airFrame.samples[i,:][0])
        
        data_i = C172.denormalize(np.array(airFrame.samples[i,:][0])).T
        try:
            LD_interpolated = f_interpolation(data_i)
        except(ValueError):
            print(data_i)
        # print(pdf, LD_interpolated)
        expected_value += LD_interpolated[0]
        pdfs.append(pdf)
    total_pdf = sum(pdfs)
    # print(total_pdf, expected_value, expected_value/len(airFrame.samples))
    expected_value = expected_value/len(airFrame.samples)
    return(expected_value)


C172 = pickle.load(open('c172_2.p', 'rb'))

data = {}
states = ['nonmorphed', 'morphed']
concepts = ['NACA0012', 'NACA4415', 'NACA641212','glider']
for state in states:
    data[state] = {}
    for concept in concepts:
        if state == 'nonmorphed':
            mat =  scipy.io.loadmat('./'+state+'/'+ state + '_' + concept)
            aoa = mat['aoa'][0]
            velocity = mat['V'][0]
            cl = mat['CL'][0]
            LD_ratio = mat['lift_to_drag']
            print(state, concept, max(LD_ratio.ravel()), expected_LD(aoa, velocity, cl, LD_ratio, C172))
        elif state == 'morphed':
            data = np.loadtxt('./'+state + '/' + state + '_' + concept + '.txt')
            aoa = data[:,0]
            velocity = data[:,1]
            cl = data[:,2]
            LD_ratio = data[:,3]
            print(state, concept, max(LD_ratio.ravel()), expected_LD(aoa, velocity, cl, LD_ratio, C172))
