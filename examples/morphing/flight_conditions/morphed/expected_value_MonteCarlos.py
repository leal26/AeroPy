import pickle
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import griddata, RegularGridInterpolator

from weather.scraper.flight_conditions import properties, Airframe


def expected(data, airFrame):
    alpha, V, lift_to_drag = data

    pdf = airFrame.pdf.score_samples(np.vstack([alpha.ravel(), V.ravel()]).T)
    pdf = np.exp(pdf.reshape(lift_to_drag.shape))
    expected_value = 0
    numerator_list = []
    denominator_list = []
    for i in range(len(lift_to_drag)):
        numerator = simps(lift_to_drag[i]*pdf[i], alpha[i])
        denominator = simps(pdf[i], alpha[i])
        numerator_list.append(numerator)
        denominator_list.append(denominator)
    numerator = simps(numerator_list, V[:, 0])
    denominator = simps(denominator_list, V[:, 0])
    expected_value = numerator/denominator
    return(expected_value)

def expected_MonteCarlos(data, airFrame):
    # data = data[data[:,0].argsort()]
    alpha = data[:,0]
    velocity = data[:,1]
    cl = data[:,2]
    lift_to_drag = data[:,3]
    expected_value = 0
    # V = 12*45
    V = 1
    pdfs = []
    f_interpolation = RegularGridInterpolator((np.unique(alpha.ravel()).T, np.unique(velocity.ravel()).T), np.reshape(lift_to_drag, [200,200]))
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
        expected_value += (V/1)*pdf[0]*LD_interpolated[0]
        pdfs.append(pdf)
    total_pdf = sum(pdfs)
    expected_value /= total_pdf
    return(expected_value)
    
C172 = pickle.load(open('c172_2.p', 'rb'))

data = {}
states = ['morphed']
concepts = ['glider', 'NACA0012', 'NACA4415', 'NACA641212']
for state in states:
    data[state] = {}
    for concept in concepts:
        if state == 'nonmorphed':
            mat =  np.loadtxt('./'+state+'/'+concept+'.txt')
            AOA, V = np.meshgrid(mat['aoa'], mat['V'])
            LD_ratio = mat['lift_to_drag']
            data[state][concept] = [AOA, V, LD_ratio]
            print(state, concept, max(LD_ratio.ravel()), expected_MonteCarlos(data[state][concept], C172)[0])
        elif state == 'morphed':
            data = np.loadtxt(state + '_' + concept + '.txt')
            LD_ratio = data[:,3]
            print(state, concept, max(LD_ratio.ravel()), expected_MonteCarlos(data, C172)[0])
