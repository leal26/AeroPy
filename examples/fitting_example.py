import aeropy.CST_3D as cst
from aeropy.filehandling.vtk import generate_surface
import aeropy.CST_3D.mesh_tools as meshtools

from scipy.spatial.distance import directed_hausdorff
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle


def convergence_study(Ny=np.linspace(10, 50, 9),
                      Ns=np.linspace(10, 50, 9)):
    Ny = Ny.astype(int)
    Ns = Ns.astype(int)
    NY, NS = np.meshgrid(Ny, Ns)
    error = np.zeros(NY.shape)
    for i in range(len(NY)):
        for j in range(len(NY[0])):
            error[i][j] = find_edges(NY[i][j], NS[i][j])
            print(i, j, error[i][j])

    f = open('error.p', 'wb')
    pickle.dump(error, f)

    fig, ax = plt.subplots()
    cs = ax.contourf(NY, NS, error)
    fig.colorbar(cs)
    plt.xlabel('Bernstein order for taper')
    plt.xlabel('Bernstein order for shear')
    plt.show()


def find_edges(Ny, Ns):
    '''
    inputs: [location, XYZ, sy, ny, xshear]'''

    XZ = [20, 10]
    A_y = (Ny+1)*[0.15]
    ny = [.5, .5]
    A_xshear = (Ns+1)*[0]

    x0 = XZ + A_y + ny + A_xshear

    solution = minimize(shape_difference, x0, args=(Ny, Ns))
    error = solution['fun']
    return error


def shape_difference(x, Ny, Ns):
    fuselage.XYZ = [x[0], Y, x[1]]
    m = 2 + Ny + 1
    fuselage.sy = cst.BernsteinPolynomial(Ny, x[2:m])
    n = m + Ns + 1
    fuselage.xshear = cst.BernsteinPolynomial(Ns, [0]+list(x[m: n]))

    # Option trying to use norm
    # x_raw, y_raw, z_raw = upper.T
    # psi_u, eta_u = fuselage.inverse(x_raw, y_raw, z_raw)
    # x_raw, y_raw, z_raw = lower.T
    # psi_l, eta_l = fuselage.inverse(x_raw, y_raw, z_raw)
    #
    # psi_e = np.vstack((np.zeros(len(eta_l)),
    #                    np.ones(len(eta_u))))
    # eta_e = np.vstack((eta_l, eta_u))
    # mesh_e = fuselage(psi_e, eta_e)
    # lower_e, upper_e = np.dstack(mesh_e)
    #
    # plt.figure()
    # plt.scatter(lower_e[:, 0], lower_e[:, 2], label='lower_e')
    # plt.scatter(upper_e[:, 0], upper_e[:, 2], label='upper_e')
    # plt.scatter(lower[:, 0], lower[:, 2], label='lower')
    # plt.scatter(upper[:, 0], upper[:, 2], label='upper')
    # plt.legend()
    # plt.show()
    # print(np.isnan(lower_e).any(), np.isnan(lower).any(), np.isnan(eta_l).any())
    # print(np.isnan(upper_e).any(), np.isnan(upper).any(), np.isnan(eta_u).any())
    # error = np.linalg.norm(lower_e - lower)
    # print(error)
    # error += np.linalg.norm(upper_e - upper)
    # print(error)

    mesh_e = fuselage(psi_e, eta_e)
    lower_e, upper_e = np.dstack(mesh_e)

    # print('XYZ', fuselage.XYZ)
    # print('sy', x[2:8])
    # print('ny', x[8:10])
    # print('xshear', [0]+list(x[10:15]))

    # plt.figure()
    # plt.scatter(lower_e[:, 0], lower_e[:, 2], label='lower_e')
    # plt.scatter(upper_e[:, 0], upper_e[:, 2], label='upper_e')
    # plt.scatter(lower[:, 0], lower[:, 2], label='lower')
    # plt.scatter(upper[:, 0], upper[:, 2], label='upper')
    # plt.legend()
    # plt.show()
    # error = np.linalg.norm(network_e-edges)

    error = directed_hausdorff(lower_e, lower)[0]
    error += directed_hausdorff(upper_e, upper)[0]
    return(error)


raw_data = pickle.load(open('fuselage.p', 'rb'))
edges = pickle.load(open('edges.p', 'rb'))
upper = pickle.load(open('upper.p', 'rb'))
lower = pickle.load(open('lower.p', 'rb'))
eta_cp = np.linspace(0, 1, 4)
# fuselage parameters

# Order
Ny = 20
Ns = 15

# Dummy start
location = [0, 0., 0]
Y = max(edges[:, 0])
XZ = [20, 10]
XYZ = [XZ[0], Y, XZ[1]]
A_y = (Ny+1)*[0.15]
A_xshear = (Ns+1)*[0]
nx = [.5, .5]
ny = [1., 0]

f_sy = cst.BernsteinPolynomial(Ny, A_y)
f_sx = cst.BernsteinPolynomial(5, [0.172802, 0.167353, 0.130747,
                                   0.172053, 0.112797, 0.168891])
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

x_raw, y_raw, z_raw = edges.T
psi_e, eta_e = meshtools.meshparameterspace((2, N_span),
                                            psi_spacing='linear',
                                            eta_spacing='cosine')
# find_edges(Ny, Ns)
convergence_study(Ny=np.linspace(10, 20, 2),
                  Ns=np.linspace(10, 20, 2))
mesh_e = fuselage(psi_e, eta_e)
network_e = np.dstack(mesh_e)
# generate mesh


psi_f, eta_f = meshtools.meshparameterspace((N_chord, N_span),
                                            psi_spacing='linear',
                                            eta_spacing='linear')

mesh_f = fuselage(psi_f, eta_f)

# format as networks
network_f = np.dstack(mesh_f)
x, y, z = network_e.T
plt.figure()
plt.scatter(x, z)
plt.scatter(x_raw, z_raw)
plt.show()
# Generating vtk files
generate_surface(network_f, "fuselage")
