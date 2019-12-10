
import matplotlib.pyplot as plt
from scipy import stats
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
g = sns.pairplot(df, diag_kind="kde")

a = df.plot.scatter(x = '$A_{l_1}$', y = '$A_{l_2}$', c='r', s=40, edgecolors=(128./255, 0, 0), alpha=0.3)
a.yaxis.label.set_size(40)
df.plot.scatter(x = '$A_{l_1}$', y = '$A_{l_3}$', c='g', s=40, edgecolors=(34./255, 85./255, 0), alpha=0.3)

df.plot.scatter(x = '$A_{l_1}$', y = '$A_{l_4}$', c='m', s=40, edgecolors=(170/255,  0, 212/255.), alpha=0.3)

df.plot.scatter(x = '$A_{l_2}$', y = '$A_{l_3}$', c='b', s=40, edgecolors=(0, 0, 128./255), alpha=0.3)

df.plot.scatter(x = '$A_{l_2}$', y = '$A_{l_4}$', c=(1, 102/255., 0), s=40, edgecolors=(170/255, 68/255., 0), alpha=0.3)

df.plot.scatter(x = '$A_{l_3}$', y = '$A_{l_4}$', c=(0, 204/255., 1), s=40, edgecolors=(0, 102/255., 128/255.), alpha=0.3)
plt.show()



print(df[(np.abs(stats.zscore(df)) > 3).all(axis=1)])
