import numpy as np
import pickle

from aeropy.filehandling.vtk import generate_points

directory = 'D:\\GitHub\\AeroPy\\examples\\JAXA_files\\processed\\'
filename = 'wing_right'
n = 100

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
f = open('wing.p', 'wb')
pickle.dump(data, f)
generate_points(data, directory+'raw_wing')

data = data[data[:, 1] == 0]
f = open('edges_w.p', 'wb')
pickle.dump(data, f)
generate_points(data, directory+'edges_w')

print(upper)
f = open('upper_w.p', 'wb')
pickle.dump(upper, f)
generate_points(upper, directory+'upper_w')

f = open('lower_w.p', 'wb')
pickle.dump(lower, f)
generate_points(lower, directory+'lower_w')
