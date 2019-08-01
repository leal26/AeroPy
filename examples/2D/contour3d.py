from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy import stats, interpolate
import seaborn as sns
import pandas as pd
import numpy as np
import pickle

data = pickle.load(open('fitting.p', 'rb'))

# list of strings
Al = np.array(data['Al'])
Al = {'$A_{l_1}$': Al[:, 1], '$A_{l_2}$': Al[:, 2],
      '$A_{l_3}$': Al[:, 3], '$A_{l_4}$': Al[:, 4]}

# Calling DataFrame constructor on list
df = pd.DataFrame(Al)
df.to_csv('coefficient_data.csv')
# g = sns.pairplot(df, diag_kind="kde")

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y, Z = axes3d.get_test_data(0.05)
ax.scatter(X, Y, Z, alpha=0.3)

N = 100
xi = np.linspace(X.min(), X.max(), N)
yi = np.linspace(Y.min(), Y.max(), N)
zi = np.linspace(Z.min(), Z.max(), N)

# Projection perpendicular to Z
values = np.array([X.flatten(), Y.flatten()])
kde = stats.gaussian_kde(values)
density = kde(values)

di = interpolate.griddata((X.flatten(), Y.flatten()), density,
                          (xi[None, :], yi[:, None]), method='cubic')
cset = ax.scatter(xi, yi, di, zdir='z', offset=-100, cmap=cm.coolwarm)

# Projection perpendicular to X
# values = np.array([Z.flatten(), Y.flatten()])
# kde = stats.gaussian_kde(values)
# density = kde(values)
#
# di = interpolate.griddata((Z.flatten(), Y.flatten()), density,
#                           (zi[None, :], yi[:, None]), method='cubic')
# print(di)
# cset = ax.contourf(zi, yi, di, zdir='x', offset=-40, cmap=cm.coolwarm)

# Projection perpendicular to Y
values = np.array([X.flatten(), Z.flatten()])
kde = stats.gaussian_kde(values)
density = kde(values)

di = interpolate.griddata((X.flatten(), Z.flatten()), density,
                          (zi[None, :], yi[:, None]), method='cubic')
print(xi.shape, di.shape, zi.shape)
cset = ax.scatter(xi, zi, di, zdir='z', offset=40, cmap=cm.coolwarm)

# plt.colorbar(label='Probability')
ax.set_xlabel('X')
ax.set_xlim(-40, 40)
ax.set_ylabel('Y')
ax.set_ylim(-40, 40)
ax.set_zlabel('Z')
ax.set_zlim(-100, 100)

plt.show()
