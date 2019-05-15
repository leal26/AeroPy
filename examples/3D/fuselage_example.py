import aeropy.CST_3D as cst
from aeropy.filehandling.vtk import generate_surface
import aeropy.CST_3D.mesh_tools as meshtools

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle

raw_data = pickle.load(open('fuselage.p', 'rb'))
edges = pickle.load(open('edges.p', 'rb'))
eta_cp = np.linspace(0, 1, 4)
# fuselage parameters
length = 1.

f_sy = cst.BernsteinPolynomial(5, [0.172802, 0.167353, 0.130747,
                                   0.172053, 0.112797, 0.168891])
f_sx = cst.BernsteinPolynomial(5, [0.172802, 0.167353, 0.130747,
                                   0.172053, 0.112797, 0.168891])
f_xshear = cst.BernsteinPolynomial(5, [0.172802, 0.167353, 0.130747,
                                       0.172053, 0.112797, 0.168891])

location = [0., 0., 0.]
nx = (.5, .5)
ny = (1., 1.)
XYZ = (.4, length, .2)

fuselage = cst.CST3D(rotation=(0., -90., -90.),
                     XYZ=XYZ,
                     nx=nx,
                     sy=f_sy,
                     ref=(0.5, 0., 0.),
                     ny=ny,
                     xshear=f_xshear)

# generate mesh
N_chord = 50
N_span = 100


psi_f, eta_f = meshtools.meshparameterspace((N_chord, N_span),
                                            psi_spacing='linear',
                                            eta_spacing='linear')

mesh_f = fuselage(psi_f, eta_f)

# format as networks
network_f = np.dstack(mesh_f)

# Generating vtk files
generate_surface(network_f, "fuselage")


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# x, y, z = network_f.T
# ax.scatter(x, y, z)
# x, y, z = raw_data.T
# ax.scatter(x, y, z)
# plt.show()
