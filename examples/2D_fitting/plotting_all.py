
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
plt.show()

print(df[(np.abs(stats.zscore(df)) > 3).all(axis=1)])
