import numpy as np
import matplotlib.pyplot as plt
import pickle
from optimization_tools.DOE import pareto_frontier

filename = 'wing_lower'
for i in range(3):
    data_i = np.genfromtxt(filename+'_%i.csv' % i, delimiter=',')
    for j in reversed(range(len(data_i))):
        if np.isnan(data_i[j]).any():
            data_i = np.delete(data_i, j, 0)
    x, y, z = data_i.T
    x_LE, y_LE, z_LE = pareto_frontier(x, y, z, False, True, tol=6e-2)
    x_TE, y_TE, z_TE = pareto_frontier(x, y, z, True, False, tol=6e-2)

    if not i:
        raw = data_i
        LE = np.vstack([x_LE, y_LE, z_LE]).T
        TE = np.vstack([x_TE, y_TE, z_TE]).T
    else:
        raw = np.vstack([raw, data_i])
        LE = np.vstack([LE, np.vstack([x_LE, y_LE, z_LE]).T])
        TE = np.vstack([TE, np.vstack([x_TE, y_TE, z_TE]).T])
plt.figure()
x, y, z = raw.T
plt.scatter(x, y)
x, y, z = LE.T
plt.scatter(x, y, c='r')
x, y, z = TE.T
plt.scatter(x, y, c='g')
plt.show()

f = open('edges', 'wb')
pickle.dump(np.array([LE, TE]), f)
f.close()

f = open('all_lower', 'wb')
pickle.dump(raw, f)
f.close()
