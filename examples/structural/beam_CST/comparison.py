import matplotlib.pyplot as plt
import numpy as np

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem


def distributed_load(s):
    try:
        len(s)
        return w*np.ones(len(s))
    except:
        return w


def format_input(input):
    # COnsidering BC for zero derivative at the root
    return list(input) + [-input[0]]


p = properties()
w = p.young*p.inertia*4/p.length**3
s = np.linspace(0, 1, 11)
l = loads(distributed_load=distributed_load, follower=False)


# Poly
poly_solution = [0., 0., -0.92054554, 0.48360499, -0.06760387]
g = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord=1, color='b')
b_poly = beam_chen(g, p, l, s)

b_poly.g.D = np.array(poly_solution)
b_poly.g.internal_variables(b_poly.length)
b_poly.g.calculate_x1(b_poly.s)
b_poly.x = b_poly.g.x1_grid
b_poly.y = b_poly.g.x3(b_poly.x)
b_poly.calculate_M()
b_poly.g.radius_curvature(b_poly.g.x1_grid)

# Cst
CST_solution = [0.475708, 0.359857, 0.250499, 0.176188]
g = CoordinateSystem.CST(D=[0, 0, 0, 0, 0], chord=1, color='b', N1=1, N2=1)
l = loads(distributed_load=distributed_load, follower=False)
b = beam_chen(g, p, l, s, ignore_ends=True)

# b.parameterized_solver(format_input=format_input, x0=[0.475708, 0.359857,
#                                                       0.250499, 0.176188])
b.g.D = np.array(CST_solution + [-CST_solution[0]])
print('D', b.g.D)
b.g.internal_variables(b.length)
b.g.calculate_x1(b.s)
b.x = b.g.x1_grid
b.y = b.g.x3(b.x)
b.calculate_M()
b.g.radius_curvature(b.g.x1_grid)
print('chord', b.g.chord)
print('deltaz', b.g.deltaz)
plt.figure()
plt.plot(b_poly.x, b_poly.M, 'b', label='poly')
plt.plot(b.x, b.M, 'r', label='CST')
plt.ylabel('Moment')
plt.legend()

plt.figure()
plt.plot(b_poly.x, b_poly.g.rho, 'b', label='poly')
plt.plot(b.x, b.g.rho, 'r', label='CST')
plt.ylabel('Radius of curvature')
plt.legend()

plt.figure()
plt.plot(b_poly.x, b_poly.y, c='b', label='poly')
plt.plot(b.x, b.y, c='r', label='CST')
plt.ylabel('y')
plt.legend()

plt.figure()
plt.scatter(b_poly.x, b_poly.g.x3(b_poly.x, 'x1'), c='b', label='poly')
plt.scatter(b.x, b.g.x3(b.x, 'x1'), c='r', label='CST')
plt.ylabel('D')
plt.legend()

plt.figure()
plt.scatter(b_poly.x, b_poly.g.x3(b_poly.x, 'x11'), c='b', label='poly')
plt.scatter(b.x, b.g.x3(b.x, 'x11'), c='r', label='CST')
plt.ylabel('DD')
plt.legend()
plt.show()
