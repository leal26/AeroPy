import aeropy.CST_3D as cst
import aeropy.CST_3D.mesh_tools as meshtools
from aeropy.filehandling.vtk import generate_surface

import time
import pickle
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from multiprocessing import Pool

raw_data = pickle.load(open('fuselage.p', 'rb'))
raw_data = raw_data[::80]
x_raw, y_raw, z_raw = raw_data.T
fuselage = pickle.load(open('fuselage_object.p', 'rb'))
print(fuselage.XYZ)
Nx = 5
n_cp = 5

N_chord = 50
N_span = 100

fuselage.nx = [0.5, 0.5]


def convergence_study(Nx=np.linspace(2, 10, 5),
                      Nc=np.linspace(2, 14, 5),
                      parallel=False):
    Nx = Nx.astype(int)
    Nc = Nc.astype(int)
    NX, NC = np.meshgrid(Nx, Nc)
    output = []

    NX_f = NX.flatten()
    NC_f = NC.flatten()

    if parallel:
        p = Pool(len(NX_f))
    else:
        p = Pool(1)

    p = Pool(len(NX_f))
    input = np.vstack([NX_f, NC_f]).T
    solutions = p.map(find_coefficients, input)
    error = np.array([solutions[i][0] for i in range(len(solutions))])
    error = error.reshape(NX.shape)

    f = open('convergence_coefficients.p', 'wb')
    pickle.dump(output, f)

    fig, ax = plt.subplots()
    cs = ax.contourf(NX, NC, error)
    fig.colorbar(cs)
    plt.xlabel('Bernstein order for shape')
    plt.ylabel('Number of control points')
    plt.show()


def find_coefficients(args, update_object=False):
    '''
    inputs: [location, XYZ, sy, ny, xshear]'''
    Nx, n_cp = args
    x0 = (1./(Nx+1))*np.ones((Nx+1)*n_cp + n_cp - 2)
    x0[:n_cp-2] = np.linspace(0, 1., n_cp)[1:-1]

    solution = minimize(shape_difference, x0, args=(Nx, n_cp))
    error = solution['fun']
    if update_object:
        _update_object(solution['x'], Nx, n_cp)
    psi_f, eta_f = meshtools.meshparameterspace((N_chord, N_span),
                                                psi_spacing='linear',
                                                eta_spacing='linear')

    mesh_f = fuselage(psi_f, eta_f)
    network_f = np.dstack(mesh_f)
    generate_surface(network_f, "fuselage_%i_%i" % (Nx, n_cp))
    print('Nx=%i, Nc=%i, error=%f' % (Nx, n_cp, error))
    return error, solution['x']


def shape_difference(x, Nx, n_cp):
    _update_object(x, Nx, n_cp)
    # Option trying to use norm
    mesh_f = convert_points()

    error = np.linalg.norm(mesh_f - raw_data)

    # error = directed_hausdorff(lower_e, lower)[0]
    # error += directed_hausdorff(upper_e, upper)[0]
    # print(error)
    return(error)


def _update_object(x, Nx, n_cp):
    eta_cp = [0] + list(x[:n_cp-2]) + [1]
    Ax_mat = x[n_cp-2:].reshape(Nx+1, n_cp)
    Ax = []
    for i in range(Nx+1):
        # Ax.append(cst.piecewise_linear(eta_cp, Ax_mat[i]))
        Ax.append(cst.BernsteinPolynomial(len(eta_cp)-1, Ax_mat[i]))
    f_sx = cst.BernsteinPolynomial(Nx, Ax)
    fuselage.sx = f_sx


def convert_points():
    psi_f, eta_f = fuselage.inverse(x_raw, y_raw, z_raw)
    psi_f[psi_f < 0] = 0
    eta_f[eta_f < 0] = 0
    psi_f[psi_f > 1] = 1
    eta_f[eta_f > 1] = 1
    # plt.figure()
    # plt.scatter(psi_f, eta_f)
    # plt.show()
    mesh_f = fuselage(psi_f, eta_f)
    return(np.dstack(mesh_f))


if __name__ == "__main__":

    convergence_study(Nx=np.linspace(10, 20, Nx),
                      Nc=np.linspace(6, 30, n_cp),
                      parallel=True)
    # error, x = find_coefficients(Nx, n_cp, update_object=True)
    # print(error)
    # f = open('fuselage_object.p', 'wb')
    # pickle.dump(fuselage, f)

    # eta_cp = [0., 0.276103, .551988, 1.]
    # x0 = np.zeros((Nx+1)*n_cp)
    # Ax_mat = x0.reshape(Nx+1, n_cp)
    # Ax = []
    # for i in range(Nx+1):
    #     Ax.append(cst.piecewise_linear(eta_cp, Ax_mat[i]))
    # f_sx = cst.BernsteinPolynomial(Nx, Ax)
    # nx = (.5, .5)
    # fuselage.sx = f_sx
    # fuselage.nx = nx

    # format as networks
    network_f = convert_points()

    x, y, z = network_f.T
    plt.figure()
    plt.scatter(x_raw, z_raw, label='raw')
    plt.scatter(x, z, label='fit')
    plt.legend()
    plt.show()
    # Generating vtk files
    psi_f, eta_f = meshtools.meshparameterspace((N_chord, N_span),
                                                psi_spacing='linear',
                                                eta_spacing='linear')

    mesh_f = fuselage(psi_f, eta_f)
    network_f = np.dstack(mesh_f)
    generate_surface(network_f, "fuselage_full")
