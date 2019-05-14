# import aeropy.CST_lib as cst
import aeropy.CST_3D as cst
import aeropy.CST_3D.mesh_tools as meshtools
import panairwrapper
from aeropy.filehandling.vtk import generate_surface
from aeropy.geometry.fitting import fitting

import pickle
import numpy as np
import matplotlib.pyplot as plt


f = open('edges', 'rb')
raw = pickle.load(f)
n1 = len(raw[0])
n2 = len(raw[1])
raw = np.vstack([raw[0], raw[1]])
f.close()


def x0(dummy1, dummy2):
    '''x0 = [taper, sweep, dihedral, twist]'''
    taper_0 = [1.0, 0.375317, .204817, .080035]
    sweep_0 = [0., 14.4476, 19.1612, 24.79405]
    dihedral_0 = [0., 0.095311, 0.172374]
    twist_0 = [0.1, -0.05, -0.1, -0.1]
    location_z_0 = [0.]
    x0 = taper_0 + sweep_0 + dihedral_0 + twist_0 + location_z_0
    return x0


def update(wing, x, Nx, n_cp):
    eta_cp = [0., 0.276004, .552007, 1.]
    eta_dihedral = [0., .552007, 1.]
    f_sy = cst.piecewise_linear(eta_cp, x[:4])
    f_xshear = cst.piecewise_linear(eta_cp, x[4:8])
    f_zshear = cst.piecewise_linear(eta_dihedral, x[8:11])
    f_twist = cst.piecewise_linear(eta_cp, x[11:15])
    location_z = x[15]

    wing.sy = f_sy
    wing.xshear = f_xshear
    wing.zshear = f_zshear
    wing.twist = f_twist
    wing.location[2] = location_z

# def calculate_points(wing, raw):
#     x_raw, y_raw, z_raw = raw.T
#     psi_w, eta_w = wing.inverse(x_raw, y_raw, z_raw)
#     psi_w[psi_w < 0] = 0
#     eta_w[eta_w < 0] = 0
#     psi_w[psi_w > 1] = 1
#     eta_w[eta_w > 1] = 1
#     mesh_w = wing(psi_w, eta_w)
#     return(np.dstack(mesh_w))


def calculate_points(wing, raw):
    x_raw, y_raw, z_raw = raw.T
    eta = wing._inverse_y(y_raw)

    psi = np.hstack((np.zeros(n1),
                     np.ones(n2)))
    eta = np.array(eta)
    psi[psi < 0] = 0
    eta[eta < 0] = 0
    psi[psi > 1] = 1
    eta[eta > 1] = 1
    mesh_e = wing(psi, eta)
    # print(np.dstack(mesh_e).shape)
    # from mpl_toolkits.mplot3d import Axes3D
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    #
    # ax.scatter(np.dstack(mesh_e)[0, :, 0], np.dstack(mesh_e)[0, :, 1],
    #            np.dstack(mesh_e)[0, :, 2])
    # ax.scatter(x_raw, y_raw, z_raw)
    # plt.show()
    return(np.dstack(mesh_e))


def generate_vtk(object, *args):
    try:
        N_chord, N_span = args
    except(ValueError):
        N_chord, N_span = [20, 50]
    psi, eta = meshtools.meshparameterspace((N_chord, N_span),
                                            psi_spacing='cosine',
                                            eta_spacing='cosine')
    print(object)
    mesh = object(psi, eta)
    network_wu = np.dstack(mesh)
    generate_surface(network_wu, "wing_upper")


if __name__ == "__main__":
    # Load raw data

    # wing_points = np.genfromtxt('./wing_upper_points.csv', skip_header=1, delimiter=',')
    aoa = 0  # 2.3067
    location_z = 0  # -0.535
    # wing parameters
    span = 4.578*2.  # meters
    eta_cp = [0., 0.276004, .552007, 1.]
    taper = [1.0, 0.375317, .204817, .080035]
    chord_root = 21.43
    sweep = [0., 14.4476, 19.1612, 24.79405]
    # sweep = [0., 0.]
    dihedral = [0., 0.095311, 0.172374]  # , 0., 0.]
    eta_dihedral = [0., .552007, 1.]
    twist = [0.1, -0.05, -0.1, -0.1]

    f_sy = cst.piecewise_linear(eta_cp, taper)

    # shape coefficients upper wing
    a_mat_upper = np.array([[0.01263563, 0.02351035, 0.01410948, 0.02870071, 0.00672103, 0.01423143],
                            [0.01493258, 0.02318954, 0.0256106, 0.02404792, 0.0188881, 0.01735003],
                            [0.01451941, 0.03214948, 0.03073803, 0.03342343, 0.03707097, 0.03925866],
                            [0.01687441, 0.01973094, 0.05100545, 0.01154137, 0.05669724, 0.02297795]]).T

    A0_u = cst.piecewise_linear(eta_cp, a_mat_upper[0])
    A1_u = cst.piecewise_linear(eta_cp, a_mat_upper[1])
    A2_u = cst.piecewise_linear(eta_cp, a_mat_upper[2])
    A3_u = cst.piecewise_linear(eta_cp, a_mat_upper[3])
    A4_u = cst.piecewise_linear(eta_cp, a_mat_upper[4])
    A5_u = cst.piecewise_linear(eta_cp, a_mat_upper[5])

    A_upper = [A0_u, A1_u, A2_u, A3_u, A4_u, A5_u]

    f_sx_upper = cst.BernsteinPolynomial(5, A_upper)

    f_xshear = cst.piecewise_linear(eta_cp, sweep)
    f_zshear = cst.piecewise_linear(eta_dihedral, dihedral)
    f_twist = cst.piecewise_linear(eta_cp, twist)

    wing_upper = cst.CST3D(rotation=(0., aoa, 0.),
                           location=[12., 0., location_z],
                           XYZ=(chord_root, span/2., chord_root),
                           ref=(0., 0., 0.),
                           sx=f_sx_upper,
                           nx=(0.5, 1.),
                           sy=f_sy,
                           ny=(0., 0.),
                           xshear=f_xshear,
                           zshear=f_zshear,
                           twist=f_twist)

    study = fitting(object=wing_upper,
                    update=update,
                    x0=x0,
                    raw=raw,
                    calculate_points=calculate_points,
                    callback=generate_vtk)
    x = study.find()['x']

    study.update(x, None, None)
    f = open('wing.p', 'wb')
    pickle.dump(study.object, f)
    study.callback(20, 50)
