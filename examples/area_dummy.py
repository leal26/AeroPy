import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.interpolate import interp2d
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


# def geometry(y, x):
#     return np.sin(y)*np.sin(x)


def mach_cone(y, MACH=1.6, y0=0):
    MACH = 1.6
    a = -1
    b = y0
    return a*y + b


def diff(y, MACH, y0, x):
    return abs(mach_cone(y, MACH, y0) - geometry(y, x))


mach = 1.6

# Create mesh
n = 10
x = np.linspace(0, np.pi, n)
y = np.linspace(0, np.pi, n)
X, Y = np.meshgrid(x, y)
Z = np.sin(X)*np.sin(Y)
geometry = interp2d(Y, X, Z, kind='cubic')
y0_list = np.linspace(0, np.pi, 50)
x_solution = np.linspace(0, np.pi)
y_solution = np.zeros(x_solution.shape)
z_solution = np.zeros(x_solution.shape)
output = []
A = []
for j in range(len(y0_list)):
    y0 = y0_list[j]
    for i in range(len(x_solution)):
        y_solution[i] = fsolve(diff, args=(mach, y0, x_solution[i]), x0=0)
        z_solution[i] = geometry(y_solution[i], x_solution[i])
    points = np.array([x_solution, y_solution, z_solution]).T
    output.append(points)
    points = points[points[:, 0].argsort()]
    # print(points)
    A.append(area(points))
    # print(area(points))
output = np.array(output)
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

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(X, Y, Z, c='b')
x, y, z = output.reshape(len(y0_list)*len(x_solution), 3).T
ax.scatter(x, y, z, c='r')
plt.show()
