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
            print('Ny=%i, Ns=%i, error=%f' % (i, j, error[i][j]))

    f = open('convergence_edges.p', 'wb')
    pickle.dump(error, f)

    fig, ax = plt.subplots()
    cs = ax.contourf(NY, NS, error)
    fig.colorbar(cs)
    plt.xlabel('Bernstein order for taper')
    plt.ylabel('Bernstein order for shear')
    plt.show()


def find_edges(Ny, Ns):
    '''
    inputs: [location, XYZ, sy, ny, xshear]'''

    XZ = [20, 10]
    A_y = (Ny+1)*[1./(Ny+1)]
    A_xshear = (Ns+1)*[1./(Ns+1)]

    x0 = XZ + A_y + A_xshear

    solution = minimize(shape_difference, x0, args=(Ny, Ns))
    error = solution['fun']
    return error


def shape_difference(x, Ny, Ns):
    fuselage.XYZ = [x[0], Y, x[1]]
    m = 2 + Ny + 1
    fuselage.sy = cst.BernsteinPolynomial(Ny, x[2:m])
    n = m + Ns

    fuselage.xshear = cst.BernsteinPolynomial(Ns, [0] + list(x[m: n]))

    # Option trying to use norm
    lower_e, upper_e = convert_points()

    error = np.linalg.norm(lower_e - lower)
    error += np.linalg.norm(upper_e - upper)
    return(error)


def convert_points():
    x_raw, y_raw, z_raw = upper.T
    psi_u, eta_u = fuselage.inverse(x_raw, y_raw, z_raw)
    x_raw, y_raw, z_raw = lower.T
    psi_l, eta_l = fuselage.inverse(x_raw, y_raw, z_raw)

    psi_e = np.vstack((np.zeros(len(eta_l)),
                       np.ones(len(eta_u))))
    eta_e = np.vstack((eta_l, eta_u))
    mesh_e = fuselage(psi_e, eta_e)
    lower_e, upper_e = np.dstack(mesh_e)
    return(lower_e, upper_e)


raw_data = pickle.load(open('fuselage.p', 'rb'))

upper = pickle.load(open('upper.p', 'rb'))
lower = pickle.load(open('lower.p', 'rb'))
eta_cp = np.linspace(0, 1, 4)
# fuselage parameters

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

find_edges(Ny, Ns)
f = open('fuselage_object.p', 'wb')
pickle.dump(fuselage, f)
network_e = convert_points()
lower_e, upper_e = network_e
# convergence_study(Ny=np.linspace(5, 25, 5),
#                   Ns=np.linspace(5, 25, 5))

plt.figure()
plt.scatter(lower_e[:, 0], lower_e[:, 2], c='b', label='fit')
plt.scatter(upper_e[:, 0], upper_e[:, 2], c='b')
plt.scatter(lower[:, 0], lower[:, 2], c='r', label='raw')
plt.scatter(upper[:, 0], upper[:, 2], c='r')
plt.legend()
plt.show()

# Generating vtk files
network_e = np.array([lower_e, upper_e])
generate_surface(network_e, "fuselage_edges")
