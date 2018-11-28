import numpy as np
import pickle

from aeropy.filehandling.vtk import generate_points

directory = 'D:\\GitHub\\AeroPy\\examples\\JAXA_files\\processed\\'
filename = 'fuselage'
n = 200

data = []
upper = []
lower = []
for i in range(n):
    print(i)
    data_i = np.genfromtxt(directory+filename+'_%i.csv' % i, delimiter=',')
    for j in reversed(range(len(data_i))):
        if np.isnan(data_i[j]).any():
            data_i = np.delete(data_i, j, 0)
    data += list(data_i)
    print(data_i[np.argmax(data_i[:, 2])])
    upper += [list(data_i[np.argmax(data_i[:, 2])],)]
    lower += [list(data_i[np.argmin(data_i[:, 2])],)]
data = np.array(data)
upper = np.array(upper)
lower = np.array(lower)

data = data[data[:, 1] >= 0]
f = open('fuselage.p', 'wb')
pickle.dump(data, f)
generate_points(data, directory+'raw_fuselage')

data = data[data[:, 1] == 0]
f = open('edges.p', 'wb')
pickle.dump(data, f)
generate_points(data, directory+'edges')

print(upper)
f = open('upper.p', 'wb')
pickle.dump(upper, f)
generate_points(upper, directory+'upper')

f = open('lower.p', 'wb')
pickle.dump(lower, f)
generate_points(lower, directory+'lower')
