import aeropy.CST_3D as cst
import aeropy.CST_3D.mesh_tools as meshtools
from aeropy.filehandling.vtk import generate_surface
from aeropy.geometry.fitting import fitting

import time
import pickle
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from multiprocessing import Pool


def x0(Nx, n_cp):
    '''
    inputs: [location, XYZ, sy, ny, xshear]'''

    x0 = (1./(Nx+1))*np.ones((Nx+1)*n_cp + n_cp - 2)
    x0[:n_cp-2] = np.linspace(0, 1., n_cp)[1:-1]
    return x0


def update(fuselage, x, Nx, n_cp):
    eta_cp = [0] + list(x[:n_cp-2]) + [1]
    Ax_mat = x[n_cp-2:].reshape(Nx+1, n_cp)
    Ax = []
    for i in range(Nx+1):
        # Ax.append(cst.piecewise_linear(eta_cp, Ax_mat[i]))
        Ax.append(cst.BernsteinPolynomial(len(eta_cp)-1, Ax_mat[i]))
    f_sx = cst.BernsteinPolynomial(Nx, Ax)
    fuselage.sx = f_sx


def calculate_points(fuselage, raw):
    x_raw, y_raw, z_raw = raw.T
    psi_f, eta_f = fuselage.inverse(x_raw, y_raw, z_raw)
    psi_f[psi_f < 0] = 0
    eta_f[eta_f < 0] = 0
    psi_f[psi_f > 1] = 1
    eta_f[eta_f > 1] = 1
    mesh_f = fuselage(psi_f, eta_f)
    return(np.dstack(mesh_f))


def generate_vtk(fuselage, Nx, n_cp, N_chord=50, N_span=100):
    psi, eta = meshtools.meshparameterspace((N_chord, N_span),
                                            psi_spacing='linear',
                                            eta_spacing='linear')

    mesh = fuselage(psi, eta)
    network = np.dstack(mesh)
    generate_surface(network, "fuselage_%i_%i" % (Nx, n_cp))


if __name__ == "__main__":
    raw = pickle.load(open('fuselage.p', 'rb'))
    raw = raw[::80]
    x_raw, y_raw, z_raw = raw.T
    fuselage = pickle.load(open('fuselage_object.p', 'rb'))
    Nx = 4
    n_cp = 4

    fuselage.nx = [0.5, 0.5]
    study = fitting(object=fuselage,
                    update=update,
                    p1=np.linspace(10, 25, Nx),
                    p2=np.linspace(10, 25, n_cp),
                    p1_name='Berstein order for Ay',
                    p2_name='Berstein order for Sx',
                    x0=x0,
                    raw=raw,
                    calculate_points=calculate_points,
                    callback=generate_vtk)
    study.convergence_study(parallel=True)
    study.plot_study()

    # Generating vtk files
    N_chord = 50
    N_span = 100
    psi_f, eta_f = meshtools.meshparameterspace((N_chord, N_span),
                                                psi_spacing='linear',
                                                eta_spacing='linear')

    mesh_f = fuselage(psi_f, eta_f)
    network_f = np.dstack(mesh_f)
    generate_surface(network_f, "fuselage_full")
