import aeropy.CST_3D as cst
from aeropy.filehandling.vtk import generate_surface
import aeropy.CST_3D.mesh_tools as meshtools
from aeropy.geometry.fitting import fitting

from scipy.spatial.distance import directed_hausdorff
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from multiprocessing import Pool
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle


def x0(Ny, Ns):
    '''
    inputs: [location, XYZ, sy, ny, xshear]'''

    XZ = [20, 10]
    A_y = (Ny+1)*[1./(Ny+1)]
    A_xshear = (Ns+1)*[1./(Ns+1)]
    x0 = XZ + A_y + A_xshear
    return x0


def update(fuselage, x, p1_i, p2_i):
    fuselage.XYZ[0] = x[0]
    fuselage.XYZ[2] = x[1]
    m = 2 + p1_i + 1
    fuselage.sy = cst.BernsteinPolynomial(p1_i, x[2: m])
    n = m + p2_i
    fuselage.xshear = cst.BernsteinPolynomial(p2_i, [0] + list(x[m: n]))


def calculate_points(fuselage, raw):
    x_raw, y_raw, z_raw = raw[0].T
    psi_u, eta_u = fuselage.inverse(x_raw, y_raw, z_raw)
    x_raw, y_raw, z_raw = raw[1].T
    psi_l, eta_l = fuselage.inverse(x_raw, y_raw, z_raw)
    psi_e = np.vstack((np.zeros(len(eta_l)),
                       np.ones(len(eta_u))))
    eta_e = np.vstack((eta_l, eta_u))
    mesh_e = fuselage(psi_e, eta_e)
    return(np.dstack(mesh_e))


if __name__ == "__main__":
    raw_data = pickle.load(open('fuselage.p', 'rb'))

    upper = pickle.load(open('upper.p', 'rb'))
    lower = pickle.load(open('lower.p', 'rb'))
    raw = np.array([lower, upper])
    eta_cp = np.linspace(0, 1, 4)

    # Order
    Ny = 30
    Ns = 20

    # Dummy start
    location = [0, 0., 0]
    Y = max([max(upper[:, 0]), max(lower[:, 0])])
    XZ = [20, 10]
    XYZ = [XZ[0], Y, XZ[1]]
    A_y = (Ny+1)*[0.15]
    A_xshear = (Ns+1)*[0]
    nx = [.5, .5]
    ny = [0.75, 0]

    f_sy = cst.BernsteinPolynomial(Ny, A_y)
    f_sx = cst.BernsteinPolynomial(1, [1, 1])
    f_xshear = cst.BernsteinPolynomial(Ns, A_xshear)

    fuselage = cst.CST3D(rotation=(0., -90., -90.),
                         location=location,
                         XYZ=XYZ,
                         nx=nx,
                         sy=f_sy,
                         ref=(0.5, 0., 0.),
                         ny=ny,
                         xshear=f_xshear)

    N_chord = 50
    N_span = 100

    # find_edges(Ny, Ns)
    # f = open('fuselage_object.p', 'wb')
    # pickle.dump(fuselage, f)

    study = fitting(object=fuselage,
                    update=update,
                    p1=np.linspace(5, 30, 6),
                    p2=np.linspace(5, 30, 6),
                    p1_name='Berstein order for Sy',
                    p2_name='Berstein order for shear',
                    x0=x0,
                    raw=raw,
                    calculate_points=calculate_points)
    study.convergence_study(parallel=True)
    study.plot_study()

    # Store convergence data
    f = open('convergence_edges.p', 'wb')
    pickle.dump(study, f)

    # Plot it
    study.plot_fit()

    # Generating vtk files
    network = calculate_points(fuselage, raw)
    generate_surface(network, "fuselage_edges")
