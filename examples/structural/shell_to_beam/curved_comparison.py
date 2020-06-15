from aeropy.geometry.parametric import CoordinateSystem
from aeropy.structural.shell import shell
from aeropy.structural.stable_solution import (mesh_1D, properties,
                                               boundary_conditions,
                                               euler_bernoulle)
from optimization_tools.DOE import DOE
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math


# bc = boundary_conditions(distributed_load=1)

# Define parent and child geometries


# Indian paper
# chord_parent = 0.4385
# width = 0.0385
# thickness = 0.00625
# experiment = {0: [0.,0.,0.5492,0.2000],
#               138.321: [0., 0., 0.3746, 0.1753],
#               211.896: [0., 0., 0.2475, 0.2864],}
#               # 287.433: [0., 0., 0.1287, 0.3582],
#               # 460.089: [0., 0., -0.1891, 0.6087],}
# bp = properties(young=200e9, dimensions=[width, thickness])
# exp_F = [0.000, 20.9634, 83.7474, 159.2844, 232.8594, 308.3964, 382.9542, 459.4704, 481.0524]
# exp_delta = [0, 0.0047, 0.0135, 0.03098, 0.0451, 0.0622, 0.0793, 0.0977, 0.103]
# bounds_opt = np.array([[-0.2,1.0],
#                    [0,1.0]])

# Chen
chord_parent = 0.4
width = 0.025
thickness = 0.0004
experiment = {0: [0, 0,-1.05540301,1.19321522],
              0.098: [0, 0, -1.72387251, 1.5985958],
              0.196: [0, 0, -2.3139154, 1.77815251],}

bp = properties(young=194.3e9, dimensions=[width, thickness])
exp_F = [0.000, 0.098, 0.196, 0.294, 0.392, 0.490, 0.588]
exp_delta = [0.089, 0.149, 0.195, 0.227, 0.251, 0.268, 0.281]

bounds_opt = np.array([[-3,-1.0],
                  [1,2]])

curve_parent = CoordinateSystem.polynomial(experiment[0], chord_parent, 'b')

colors = ['0.0','0.3', '0.5', '0.7']
loads = list(experiment.keys())
span_num = []

for i in range(len(loads)):
    load = loads[i]
    bc = boundary_conditions(concentrated_load=np.array([[0, 0, -load], ]), load_x = [chord_parent])
    curve_child = CoordinateSystem.polynomial(experiment[load], chord_parent, color = 'g')

    # Define shell
    beam = shell(curve_parent, curve_child, bp, bc)

    # Kinematics
    bounds = np.array([[0, 1.5*chord_parent],])
    beam.calculate_chord(bounds = bounds)
    beam.theta1 = np.linspace(0, beam.g_p.arclength()[0], 20)
    beam.g_p.bounds = bounds
    beam.g_c.bounds = bounds

    beam.update_parent()
    beam.update_child()
    plt.figure('Results')
    if i == 0:
        beam.g_p.plot(label='Parent', linestyle = '-')
    else:
        beam.g_c.plot(label='Experiment: %f N' % load, linestyle = '-',
                      color = colors[i-1], scatter = True)

        beam.minimum_potential(x0=experiment[0][2:],
                               input_function = lambda x: [0,0] + list(x),
                               bounds = bounds_opt)
        # coefficients, results = beam.stepped_loading(x0=experiment[0][2:],
        #                                              input_function = lambda x: [0,0] + list(x),
        #                                              bounds = np.array([[-0.2,1.0],
        #                                                                 [0,1.0]]))
                                                     # bounds = np.array([[-3,-1],
                                                     #                    [1,2]]))


        beam.g_c.plot(label='Minimum: %f N' % load, linestyle = '-', color = colors[i-1])

        span_num.append(beam.g_c.chord)

    # plt.figure('summary')
    # plt.plot(-1*np.array(results[:,0]), -1*np.array(results[:,1]), colors[i-1])

    # x, R = beam.minimum_potential(x0 = beam.g_p.D, bounds = [[.5*beam.g_p.D, 2*beam.g_p.D],])
    # print(x, R)
    # plt.figure()
    # beam.g_p.plot(label='Parent')
    # # beam.g_c.plot(label='Euler-Bernoulle', color = 'k')
    # # beam.g_c.D[2:] = x
    # beam.update_child()
    # beam.g_c.plot(label='Minimum Potential', color = '.5', linestyle= '--')
    # plt.legend()


# plt.figure('summary')
# plt.scatter(exp_F, exp_delta, c='g', label='Experiments', edgecolors='k', zorder = 10)
plt.legend()
plt.show()
BRAKE
def DOE_function(inputs):
    beam.g_c.D[2] = inputs['D2']
    beam.g_c.D[3] = inputs['D3']
    beam.update_child()
    return({'U': beam.U, 'W':beam.W, 'R':beam.R})


# Define points
print(eulerBernoulle)
problem = DOE(levels=20, driver='Full Factorial')
problem.add_variable('D2', lower=0, upper=2*eulerBernoulle.g.D[2], type=float)
problem.add_variable('D3', lower=2*eulerBernoulle.g.D[3], upper=0, type=float)
problem.define_points()

# Run for a function with dictionary as inputs
problem.run(DOE_function, tracking=True)

# Plot factor effects
# problem.plot(xlabel=['D2', 'D3'],
#              ylabel=['Strain Energy'])
# Plot domain
eulerBernoulle.analytical_solutions()
print(eulerBernoulle.g.D)
problem.plot_contour('D2', 'D3', 'W',  labels=['$D_2$', '$D_3$', 'Work'])
plt.scatter(eulerBernoulle.g.D[2], eulerBernoulle.g.D[3], marker = '^',
            color='white', label='Euler-bernoulli')
plt.scatter(x[0], x[1], marker = 's', color='k', label='Minimum Potential')
plt.legend()

problem.plot_contour('D2', 'D3', 'U',  labels=['$D_2$', '$D_3$', 'Strain Energy'])
plt.scatter(eulerBernoulle.g.D[2], eulerBernoulle.g.D[3], marker = '^',
            color='white', label='Euler-bernoulli')
plt.scatter(x[0], x[1], marker = 's', color='k', label='Minimum Potential')
plt.legend()

problem.plot_contour('D2', 'D3', 'R',  labels=['$D_2$', '$D_3$', 'Residual'])
plt.scatter(eulerBernoulle.g.D[2], eulerBernoulle.g.D[3], marker = '^',
            color='white', label='Euler-bernoulli')
plt.scatter(x[0], x[1], marker = 's', color='k', label='Minimum Potential')
plt.legend()
plt.show()
