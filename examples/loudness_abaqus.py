import numpy as np
import pickle
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, minimize, differential_evolution
from scipy.interpolate import interp2d, CloughTocher2DInterpolator, LinearNDInterpolator
from scipy.spatial import Delaunay

from rapidboom import AxieBump
from weather.boom import read_input
from weather.scraper.twister import process_data
import platform

'''Works for:
 - any cross section calculation for any angle
 - from data interpolate a cubic spline'''


def area(pts):
    'Area of cross-section.'

    if list(pts[0]) != list(pts[-1]):
        pts = pts + pts[:1]
    # x = pts[:, 0]
    # y = pts[:, 1]
    # z = pts[:, 2]
    pts[:, 1] = pts[:, 1] - pts[0, 1]
    s = 0
    for i in range(len(pts) - 1):
        s -= np.cross(pts[i, :], pts[i+1, :])
    return np.linalg.norm(s)/2


def calculating_area(X, Y, Z, y0_list, nx):
    def diff(y, MACH, y0, x):
        return abs(mach_cone(y, MACH, y0) - geometry([y, x]))

    geometry = LinearNDInterpolator(
        np.array([Y.ravel(), X.ravel()]).T, Z.ravel(), rescale=True)
    # def geometry(xi):
    #     return griddata(np.array([Y.ravel(), X.ravel()]).T, Z.ravel(), xi, method='cubic')
    x_solution = np.linspace(x_min, x_max, nx)
    y_solution = np.zeros(x_solution.shape)
    z_solution = np.zeros(x_solution.shape)
    output = []
    A = []
    for j in range(len(y0_list)):
        y0 = y0_list[j]
        for i in range(len(x_solution)):
            y_solution[i] = fsolve(diff, args=(mach, y0, x_solution[i]), x0=y0)
            z_solution[i] = geometry([y_solution[i], x_solution[i]])
        points = np.array([x_solution, y_solution, z_solution]).T
        output.append(points)
        points = points[points[:, 0].argsort()]
        A.append(area(points))

    return np.array(A), np.array(output)


def calculate_radius(y, X, Y, Z, nx, A0):
    A, output = calculating_area(X, Y, Z, y, nx)
    all_output.append(output)
    sign = np.sign(A-A0)
    r = np.sqrt(sign*(A-A0)/np.pi)
    print('A', A-A0)
    plt.plot(y, A)
    # print('r', r)
    return sign*r


def mach_cone(y, MACH, y0):
    a = -1/MACH
    b = -a*y0
    return a*y + b


def calculate_loudness(bump_function):
    # Bump design variables
    location = 12.5
    width = 2
    # Flight conditions inputs
    alt_ft = 50000.

    # Setting up for
    CASE_DIR = "./"  # axie bump case
    PANAIR_EXE = 'panair.exe'
    SBOOM_EXE = 'sboom_windows.dat.allow'

    # Run
    # axiebump = AxieBump(CASE_DIR, PANAIR_EXE, SBOOM_EXE) # for standard atmosphere
    axiebump = AxieBump(CASE_DIR, PANAIR_EXE, SBOOM_EXE, altitude=alt_ft,
                        deformation='custom')
    axiebump.MESH_COARSEN_TOL = 0.00045
    axiebump.N_TANGENTIAL = 20
    loudness = axiebump.run([bump_function, location, width])

    return loudness


all_output = []
mach = 1.6
nx = 50
ny = 20
f = open('outputs_simple_dimple.p', 'rb')  #
data = pickle.load(f, encoding='latin1')

Z, X, Y = data['COORD']['Step-2'][0].T
U3, U1, U2 = data['U']['Step-2'][0].T
Z = -Z
U3 = -U3

x_min, x_max = min(X), max(X)
y_min, y_max = min(Y), max(Y)
z_min, z_max = min(Z), max(Z)

# calculate original area
dY = (max(Y) - min(Y))
X0 = np.concatenate((X[:-1], X[:-1], X + U1, X[1:]))
Y0 = np.concatenate((Y[:-1] - 2*dY + .5, Y[:-1] - dY + .5, Y + .5 + U2, Y[1:] + dY + .5))
Z0 = np.concatenate((Z[:-1], Z[:-1], Z + U3, Z[1:]))
A0, output0 = calculating_area(X0, Y0, Z0, [min(Y)], nx)
print(A0)
A0 = A0[0]
steps = ['Step-2', 'Step-3']
loudness = {}
plt.figure()
for step in steps:
    loudness[step] = []
    for i in range(len(data['COORD'][step])):
        Z, X, Y = data['COORD'][step][i].T
        U3, U1, U2 = data['U'][step][i].T
        Z = -Z
        U3 = - U3
        # Calculate morphed area
        X = np.concatenate((X[:-1], X[:-1], X+U1, X[1:]))
        Y = np.concatenate((Y[:-1] - 2*dY + .5, Y[:-1] - dY + .5, Y + U2 + .5, Y[1:] + dY + .5))
        Z = np.concatenate((Z[:-1], Z[:-1], Z + U3, Z[1:]))

        loudness_i = calculate_loudness(lambda xx: calculate_radius(xx-12.5,
                                                                    X=X, Y=Y,
                                                                    Z=Z, nx=nx,
                                                                    A0=A0))
        loudness[step].append(loudness_i)
        print(step, i, loudness_i)
f = open('loudness.p', 'wb')
pickle.dump(loudness, f)
f.close()
f = open('output.p', 'wb')
pickle.dump(all_output, f)
f.close()
plt.show()
plt.figure()
plt.plot(data['Time']['Step-2'][:len(loudness['Step-2'])],
         loudness['Step-2'], label='Heating')
plt.plot(data['Time']['Step-3'][:len(loudness['Step-3'])],
         loudness['Step-3'], label='Cooling')
plt.legend()
plt.show()


y0_list = np.linspace(-1.5, 2, ny)
A, output = calculating_area(X, Y, Z, y0_list, nx)
A = A - A0

plt.figure()
plt.plot(y0_list, A)
plt.ylabel('Area along Mach Cone')
plt.xlabel('Distance along aircraft')
plt.show()

fig = plt.figure()
ax = Axes3D(fig)
# ax.scatter(X, Y, Z, c='b')
x, y, z = output.reshape(nx*ny, 3).T
ax.scatter(x, y, z, c='r')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
