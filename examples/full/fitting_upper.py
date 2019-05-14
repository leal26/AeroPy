# import aeropy.CST_lib as cst
import aeropy.CST_3D as cst
import aeropy.CST_3D.mesh_tools as meshtools
import panairwrapper
from aeropy.filehandling.vtk import generate_surface
from aeropy.geometry.fitting import fitting

import pickle
import numpy as np
import matplotlib.pyplot as plt


def x0_upper(dummy1, dummy2):
    '''x0 = [taper, sweep, dihedral, twist]'''
    a_mat_upper = np.array([[0.01263563, 0.02351035, 0.01410948, 0.02870071,
                             0.00672103, 0.01423143],
                            [0.01493258, 0.02318954, 0.0256106, 0.02404792,
                             0.0188881, 0.01735003],
                            [0.01451941, 0.03214948, 0.03073803, 0.03342343,
                             0.03707097, 0.03925866],
                            [0.01687441, 0.01973094, 0.05100545, 0.01154137,
                             0.05669724, 0.02297795]]).T

    x0 = a_mat_upper.flatten()
    return x0


def x0_lower(dummy1, dummy2):
    '''x0 = [taper, sweep, dihedral, twist]'''
    a_mat_lower = np.array([[-0.01115817, -0.00544143, -0.0136072, -0.01172747,
                             -0.01432724, -0.01362163],
                            [-0.01205803, 0.00105121, -0.01297246, -
                                0.00643726, -0.0084237, -0.00878171],
                            [-0.0115331, 0.00720371, 0.00013933, -0.0095523,
                             0.02167392, 0.00636361],
                            [-0.00905945, 0.0107793, -0.01631087, 0.0092086,
                             0.01001222, 0.01334978]]).T

    x0 = a_mat_lower.flatten()
    return x0


def update(wing, x, Nx, n_cp):
    eta_cp = [0., 0.276004, .552007, 1.]
    a_mat_upper = x.reshape((6, 4))

    A0_u = cst.piecewise_linear(eta_cp, a_mat_upper[0])
    A1_u = cst.piecewise_linear(eta_cp, a_mat_upper[1])
    A2_u = cst.piecewise_linear(eta_cp, a_mat_upper[2])
    A3_u = cst.piecewise_linear(eta_cp, a_mat_upper[3])
    A4_u = cst.piecewise_linear(eta_cp, a_mat_upper[4])
    A5_u = cst.piecewise_linear(eta_cp, a_mat_upper[5])
    A_upper = [A0_u, A1_u, A2_u, A3_u, A4_u, A5_u]

    wing.sx = cst.BernsteinPolynomial(5, A_upper)


def calculate_points(wing, raw):
    x_raw, y_raw, z_raw = raw.T
    psi_f, eta_f = wing.inverse(x_raw, y_raw, z_raw)
    psi_f[psi_f < 0] = 0
    eta_f[eta_f < 0] = 0
    psi_f[psi_f > 1] = 1
    eta_f[eta_f > 1] = 1
    mesh_f = wing(psi_f, eta_f)
    return(np.dstack(mesh_f))


def generate_vtk_upper(object, *args):
    try:
        N_chord, N_span = args
    except(ValueError):
        N_chord, N_span = [20, 50]
    psi, eta = meshtools.meshparameterspace((N_chord, N_span),
                                            psi_spacing='cosine',
                                            eta_spacing='cosine')

    mesh = object(psi, eta)
    network_wu = np.dstack(mesh)
    generate_surface(network_wu, "wing_upper")


def generate_vtk_lower(object, *args):
    try:
        N_chord, N_span = args
    except(ValueError):
        N_chord, N_span = [20, 50]
    psi, eta = meshtools.meshparameterspace((N_chord, N_span),
                                            psi_spacing='cosine',
                                            eta_spacing='cosine')

    mesh = object(psi, eta)
    network_wu = np.dstack(mesh)
    generate_surface(network_wu, "wing_lower")


if __name__ == "__main__":
    # Load raw data
    f = open('wing.p', 'rb')
    wing_upper = pickle.load(f)
    f.close()
    f = open('wing.p', 'rb')
    wing_lower = pickle.load(f)
    f.close()

    f = open('all_upper', 'rb')
    raw = pickle.load(f)
    f.close()

    study = fitting(object=wing_upper,
                    update=update,
                    x0=x0_upper,
                    raw=raw,
                    calculate_points=calculate_points,
                    callback=generate_vtk_upper)
    x = study.find()['x']

    f = open('wing_upper.p', 'wb')
    study.update(wing_upper, x, None, None)
    pickle.dump(wing_upper, f)
    study.callback(50, 100)

    f = open('all_lower', 'rb')
    raw = pickle.load(f)
    f.close()
    study = fitting(object=wing_lower,
                    update=update,
                    x0=x0_lower,
                    raw=raw,
                    calculate_points=calculate_points,
                    callback=generate_vtk_lower)
    x = study.find()['x']

    f = open('wing_lower.p', 'wb')
    study.update(wing_lower, x, None, None)
    pickle.dump(wing_lower, f)
    study.callback(50, 100)
