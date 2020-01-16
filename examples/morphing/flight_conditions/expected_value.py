import pickle
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

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

C172 = pickle.load(open('c172.p', 'rb'))

data = {}
states = ['nonmorphed', 'morphed']
concepts = ['UAG8814320', 'NACA0012', 'NACA4415', 'NACA641212']
for state in ['nonmorphed', 'morphed']:
    data[state] = {}
    for concept in concepts:
        if state == 'nonmorphed':
            mat =  scipy.io.loadmat('./'+state+'/'+concept+'.mat')
            AOA, V = np.meshgrid(mat['aoa'], mat['V'])
            data[state][concept] = [AOA, V, mat['lift_to_drag']]
            print(state, concept, expected(data[state][concept], C172))
        elif state == 'morphed':
            if concept == 'UAG8814320':
                mat =  scipy.io.loadmat('./'+state+'/comparisons.mat')['glider'][0][0]
            else:
                mat =  scipy.io.loadmat('./'+state+'/comparisons.mat')[concept][0][0]
            AOA = mat[-8].T
            V = mat[-7].T
            LD_ratio = mat[-2].T
            # AOA, V = np.meshgrid(mat['aoa'], mat['V'])
            data[state][concept] = [AOA, V, LD_ratio]
            # print(state, concept)
            # print(data[state][concept])
            # print(np.shape(mat['aoa']), np.shape(mat['V']), np.shape(mat['lift_to_drag']))
            print(state, concept, expected(data[state][concept], C172))
        print(AOA)

BREAK
C172 = pickle.load(open('c172.p', 'rb'))

alpha, V, lift_to_drag = owl

plt.show()
print('owl', expected(owl, C172))
print('NACA0012', expected(naca0012, C172))
print('NACA4415', expected(naca4415, C172))
