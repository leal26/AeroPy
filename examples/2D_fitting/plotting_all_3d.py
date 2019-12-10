
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from sklearn.neighbors.kde import KernelDensity
from scipy import stats
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

# Option 2: pdf

# intializing figure
fig = plt.figure(figsize=(7,5))
ax = Axes3D(fig)

# intitializing variables
points_list = [[Al['$A_{l_1}$'], Al['$A_{l_2}$']], [Al['$A_{l_1}$'], Al['$A_{l_3}$']],
               [Al['$A_{l_2}$'], Al['$A_{l_3}$']]]
direction = ['z', 'y', 'x']

# looping around axes
for i in range(len(points_list)):

    xmin = points_list[i][0].min()
    xmax = points_list[i][0].max()
    ymin = points_list[i][1].min()
    ymax = points_list[i][1].max()
    
    nx = 400
    ny = 400
    j = 1j
    # Creating pdf of dataset using scipy.stats
    One, Two = np.mgrid[-1:1:nx*j, -1:1:ny*j]  # FIXME: make mesh more fine
    positions = np.vstack([One.ravel(), Two.ravel()])
    values = np.vstack(points_list[i])
    kernel = stats.gaussian_kde(values)
    a = kernel(positions)
    fn = 0.5
    filtered = a[a>fn]
    print(i, sum(filtered)/sum(a))
    a = 1*(a>1)
    # kernel(positions).T
    a[a<fn] = 0
    a[a>=fn] = 256
    G = np.reshape(a.T, One.shape)  # FIXME: dark colors
    

    # setting directional variables
    if direction[i] == 'z':
        X = One
        Y = Two
        Z = np.ones((nx,ny)) * -1
    elif direction[i] == 'y':
        X = One
        Y = np.ones((nx,ny)) * 1
        Z = Two
    elif direction[i] == 'x':
        X = np.ones((nx,ny)) * -1
        Y = One
        Z = Two
    # Plotting contour surface on the axes
    surf = ax.plot_surface(X,Y,Z, facecolors=cm.Greys(G), alpha=1, shade=False)
    # m = cm.ScalarMappable(cmap=cm.Blues, norm=surf.norm) 
    # m.set_array(G)
    # plt.colorbar(m)  
''' # I can't get this to work
m = cm.ScalarMappable(cmap=cm.jet)
m.set_array(np.reshape(np.array(grid_u), (ny,nx)))
fig.colorbar(m)
'''
'''
# intializing figure
fig = plt.figure()
ax = Axes3D(fig)

# FIXME
# Creating pdf of dataset using sklearn
pdf = KernelDensity(kernel='gaussian').fit(values)
X, Y = np.mgrid[xmin:xmax:20j, ymin:ymax:20j]
Z = np.exp(pdf.score_samples(positions).T)
Z = np.reshape(Z, X.shape)

# Plotting contour surface on the axes
G = np.ones((20,20))
G = G * -0.75
surf = ax.plot_surface(X,Y,G, facecolors=cm.viridis(Z))
'''

ax.grid(False)
# ax.set_xticks([])
# ax.set_yticks([])
# ax.set_zticks([])

# make the panes transparent
# ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# Get rid of the spines                         
# ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
# ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
# ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

ax.set_xlim([-1.0, 1.0])
ax.set_ylim([-1.0, 1.0])
ax.set_zlim([-1.0, 1.0])

ax.set_xlabel('$A_{l_1}$')
ax.set_ylabel('$A_{l_2}$')
ax.set_zlabel('$A_{l_3}$')


fig = plt.figure(figsize=(7,5))
ax = Axes3D(fig)
# Plotting 3D scatter of points
ax.scatter(Al['$A_{l_1}$'], Al['$A_{l_2}$'], Al['$A_{l_3}$'], c='k', alpha=0.1, linewidths = 0)


ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

# make the panes transparent
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# Get rid of the spines                         
# ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
# ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
# ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

ax.set_xlim([-1.0, 1.0])
ax.set_ylim([-1.0, 1.0])
ax.set_zlim([-1.0, 1.0])

ax.set_xlabel('$A_{l_1}$')
ax.set_ylabel('$A_{l_2}$')
ax.set_zlabel('$A_{l_3}$')

plt.show()
