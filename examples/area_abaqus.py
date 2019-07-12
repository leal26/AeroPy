import numpy as np
import pickle
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.interpolate import interp2d

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
        return abs(mach_cone(y, MACH, y0) - geometry(y, x))
    geometry = interp2d(Y, X, Z, kind='cubic')
    x_solution = np.linspace(x_min, x_max, nx)
    y_solution = np.zeros(x_solution.shape)
    z_solution = np.zeros(x_solution.shape)
    output = []
    A = []
    for j in range(len(y0_list)):
        y0 = y0_list[j]
        for i in range(len(x_solution)):
            y_solution[i] = fsolve(diff, args=(mach, y0, x_solution[i]), x0=y0)
            z_solution[i] = geometry(y_solution[i], x_solution[i])
        points = np.array([x_solution, y_solution, z_solution]).T
        output.append(points)
        points = points[points[:, 0].argsort()]
        A.append(area(points))

    return np.array(A), np.array(output)


def calculate_radius(y, X, Y, Z, nx, A0):
    A, output = calculating_area(X, Y, Z, y, nx)
    return np.sqrt((A-A0)/np.pi)


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
                        deformation='Custom')
    axiebump.MESH_COARSEN_TOL = 0.00045
    axiebump.N_TANGENTIAL = 20
    loudness = axiebump.run([bump_function, location, width])

    return loudness


mach = 1.6
nx = 50
ny = 20
f = open('outputs_simple_dimple.p', 'rb')  #
data = pickle.load(f, encoding='latin1')

Z, X, Y = data['COORD']['Step-3'][0].T
U3, U1, U2 = data['U']['Step-3'][-1].T
# print(min(U3), min(U1), mi(U2))
Z = -Z
U3 = - U3
x_min, x_max = min(X), max(X)
y_min, y_max = min(Y), max(Y)
z_min, z_max = min(Z), max(Z)
# dY = (max(Y) - min(Y))
# X0 = np.concatenate((X + U1, X, X))
# Y0 = np.concatenate((Y + U2, Y + dY, Y + 2*dY))
# Z0 = np.concatenate((Z + U3, Z, Z))
# calculate original area
A0, output0 = calculating_area(X, Y, Z, [y_min], nx)
A0 = A0[0]
# Calculate morphed area
U3, U1, U2 = data['U']['Step-3'][0].T
dY = (max(Y) - min(Y))
X = np.concatenate((X + U1, X, X))
Y = np.concatenate((Y + U2, Y + dY, Y + 2*dY))
Z = np.concatenate((Z - U3, Z, Z))

y0_list = np.linspace(y_min, y_max*1.5, ny)
A, output = calculating_area(X, Y, Z, y0_list, nx)
A = A - A0

# r = calculate_radius(y0_list, X, Y, Z, nx, A0)
# y = np.linspace(0, np.pi, 100)
# m = mach_cone(y, mach, y0)
# g = geometry(y, x_solution)
#
# plt.plot(y, m, 'm')
# plt.plot(y, g, 'g')
# plt.scatter(y_solution, z_solution, c='r')
# plt.show()


# Xr = X.ravel()
# Yr = Y.ravel()
# Zr = Z.ravel()
# A = []
# for i in range(n):
#     points = np.array([Xr[Yr == y[i]], Zr[Yr == y[i]]]).T
#     points = points[points[:, 0].argsort()]
#     A.append(area(points))

plt.figure()
plt.plot(y0_list, A)
plt.ylabel('Area along Mach Cone')
plt.xlabel('Distance along aircraft')
plt.show()

# plt.figure()
# plt.plot(y0_list, r)
# plt.ylabel('Radius along Mach Cone')
# plt.xlabel('Distance along aircraft')
# plt.show()

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(X, Y, Z, c='b')
x, y, z = output.reshape(nx*ny, 3).T
ax.scatter(x, y, z, c='r')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
