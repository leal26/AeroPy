import numpy as np
import pickle
import math
import matplotlib.pyplot as plt

from airfoil_module import CST
from xfoil_module import output_reader

data = pickle.load(open('shape_study.p','rb'))

n=5
print 'For n=', n, ': '
print '\t error: ', data['error'][n-1]
print '\t deltaz: ', data['deltaz'][n-1]
print '\t Au: ', data['Au'][n-1]
print '\t Al: ', data['Al'][n-1]

filename = 'sampled_airfoil_data.csv'
raw_data = output_reader(filename, separator = ', ', header = ['x', 'y'])

# Rotating airfoil 
x_TE = (raw_data['x'][0] + raw_data['x'][-1])/2.
y_TE = (raw_data['y'][0] + raw_data['y'][-1])/2.

theta_TE = math.atan(-y_TE/x_TE)

# position trailing edge at the x-axis
processed_data = {'x':[], 'y':[]}
for i in range(len(raw_data['x'])):
    x = raw_data['x'][i]
    y = raw_data['y'][i]
    c_theta = math.cos(theta_TE)
    s_theta = math.sin(theta_TE)
    x_rotated = c_theta*x - s_theta*y
    y_rotated = s_theta*x + c_theta*y
    processed_data['x'].append(x_rotated)
    processed_data['y'].append(y_rotated)
raw_data = processed_data

# determine what is the leading edge and the rotation angle beta
processed_data = {'x':[], 'y':[]}
min_x_list = []
min_y_list = []

min_x = min(raw_data['x'])
min_index = raw_data['x'].index(min_x)
min_y = raw_data['y'][min_index]

chord = max(raw_data['x']) - min(raw_data['x'])
beta = math.atan((y_TE - min_y)/(x_TE - min_x))

for i in range(len(raw_data['x'])):
    processed_data['x'].append((raw_data['x'][i] - min_x)/chord)
    processed_data['y'].append(raw_data['y'][i]/chord)    
raw_data = processed_data

psi = np.linspace(0,1,200)
xi = CST(psi, 1., [data['deltaz'][n-1]/2., data['deltaz'][n-1]/2.], Au = data['Au'][n-1], Al= data['Al'][n-1])
plt.figure()
plt.plot(psi, xi['u'], psi, xi['l'])
plt.scatter(raw_data['x'], raw_data['y'])
n = 8
plt.xlim(0,1)
x = np.linspace(2,2*n,n)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

plt.figure()
plt.plot(x, data['error'])
plt.scatter(x, data['error'])
plt.xlabel('Number of shape functions')
plt.ylabel('Hausdorff distance (adimensional)')
plt.grid()
plt.show()