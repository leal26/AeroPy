import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import stl
from scipy.optimize import minimize, differential_evolution
from scipy.spatial.distance import directed_hausdorff
import pickle
from numpy import linalg as LA

from aeropy.xfoil_module import output_reader
from aeropy.geometry.wing import Rotation_euler_vector
from aeropy.CST_3D.module import CST_3D, dimensionalize, ControlPoints


def process_data(upper, lower, LE, TE):
    '''Displace and normalize all data and LE and TE. Only necessary if not
       processed at filehandling/stl_processing'''
    # For x
    norm = max(np.absolute(lower['y']))
    min_x = min(lower['x'])
    upper['x'] = (upper['x'] - min(upper['x']))/norm
    lower['x'] = (lower['x'] - min_x)/norm
    LE['x'] = (LE['x'] - min_x)/norm
    TE['x'] = (TE['x'] - min_x)/norm

    # For y
    upper['y'] = -1*np.array(upper['y'])/norm
    lower['y'] = -1*np.array(lower['y'])/norm
    LE['y'] = - np.array(LE['y'])/norm
    TE['y'] = - np.array(TE['y'])/norm

    # For z
    index_min = np.where(upper['x'] == min(upper['x']))[0][0]
    min_z = upper['z'][index_min]
    upper['z'] = (upper['z'] - upper['z'][index_min])/norm
    LE['z'] = (LE['z'] - min_z)/norm
    TE['z'] = (TE['z'] - min_z)/norm

    # Sort x, y, and z sort data
    [y_LE, x_LE, z_LE] = zip(*sorted(zip(LE['y'], LE['x'], LE['z'])))
    [y_TE, x_TE, z_TE] = zip(*sorted(zip(TE['y'], TE['x'], TE['z'])))

    # Correction to get 1.00000001= 1.0
    y_LE = list(y_LE)
    y_LE[-1] = 1.0
    y_LE = np.array(y_LE)
    y_TE = list(y_TE)
    y_TE[-1] = 1.0
    y_TE = np.array(y_TE)

    Raw_LE = {'x': x_LE, 'y': y_LE, 'z': z_LE}
    Raw_TE = {'x': x_TE, 'y': y_TE, 'z': z_TE}

    return upper, lower, Raw_LE, Raw_TE


def find_ControlPoints(cp, Raw_LE, Raw_TE,
                       solver='gradient'):
    '''
    inputs: [twist,  chord, sweep, shear, x_twist]'''
    m = len(cp.eta)  # number of cross sections

    # bounds: [twist,  chord, sweep, shear, x_twist]
    bounds = m*[[-1., 1.]] + m*[[.1, 4.], ] + m*[[0., 6.]] +  \
        m*[[-1., 1.]] + [[0., 4.]]

    if solver == 'differential_evolution':
        result = differential_evolution(shape_difference, bounds,
                                        disp=True, popsize=40,
                                        args=[cp, Raw_LE, Raw_TE])
        x = result.x
        error = result.fun

    if solver == 'gradient':
        x0 = m*[0., ] + m*[0., ] + m*[0.] + m*[0.] + [0.]

        solution = minimize(shape_difference, x0, bounds=bounds,
                            options={'maxfun': 30000, 'eps': 1e-07},
                            args=(cp, Raw_LE, Raw_TE))
        # print(solution)
        x = solution['x']
        error = solution['fun']

    cp = _process_input(x, cp)
    return cp, error


def shape_difference(x, cp, Raw_LE, Raw_TE, separate_errors=False):
    '''Calculates the least square error for TE AND LE.
       - x: [twist,  chord, sweep, shear, x_twist]'''

    cp = _process_input(x, cp)

    # Calculate edge coordinates
    [Data_LE, Data_TE] = edges(cp, {'LE': Raw_LE[:, 1], 'TE': Raw_TE[:, 1]})

    # Calculate error for leading edge
    error = LA.norm(Data_LE-Raw_LE)
    error += LA.norm(Data_TE-Raw_TE)

    return error


def edges(cp, eta=None):
    '''Function to model edge for trailing and leading edges'''
    if eta is None:
        N = 1000
        eta = {'LE': np.linspace(0, 1, N),
               'TE': np.linspace(0, 1, N)}
    else:
        eta = {'LE': np.array(eta['LE']),
               'TE': np.array(eta['TE'])}

    psi = {'LE': np.zeros(len(eta['LE'])),
           'TE': np.ones(len(eta['TE']))}
    xi = {'LE': np.zeros(len(eta['LE'])),
          'TE': np.zeros(len(eta['TE']))}

    Data = {}
    for key in psi:
        # Merge data
        raw_data = np.array([psi[key], eta[key], xi[key]]).T
        # Determine control points
        cp._set_functions()

        # Find solution in physical domain
        Data[key] = dimensionalize(raw_data, cp=cp, inverse=False)

    return Data['LE'], Data['TE']


def _process_input(x, cp=None, m=None, n=None, case='control_points'):
    '''Process input to and from optimizer. Used for find_ControlPoints and
       find_coefficients.

       - cp: ControlPoints object that will be updated for new values. Only
             relevant when determining ControlPoints.
       - m and n: dimensions of B coefficient matrix. Only relevant for
            case == 'coefficients'
       - case: 'control_points' and 'coefficients' '''
    # Round to a meaningful precision
    for i in range(len(x)):
        x[i] = round(x[i], 7)
    if case == 'control_points':
        m = len(cp.eta)
        cp.twist = x[:m]
        cp.chord = x[m:2*m]
        cp.sweep = x[2*m:3*m]
        cp.shear = x[3*m:4*m]
        cp.twist_origin['x'] = x[-1]
        return cp
    elif case == 'coefficients':
        nn = n+1
        mm = m+1
        Bu = np.array(x[:nn*mm]).reshape(mm, nn)
        Bl = np.array(x[nn*mm:]).reshape(mm, nn)
        return Bu, Bl


def find_coefficients(n, m, original_data, solver='gradient'):
    '''inputs: B in vector form'''
    nn = n+1  # order of Sx
    mm = m+1  # order of Sy

    # bounds: [twist,  chord, sweep, shear, x_twist]
    bounds = (nn*mm)*[[0, 1.]] + (nn*mm)*[[-1, .5]]

    if solver == 'differential_evolution':
        result = differential_evolution(shape_difference_nondimensional,
                                        bounds, disp=True, popsize=40,
                                        args=[original_data])
        x = result.x
        error = result.fun

    if solver == 'gradient':
        x0 = np.array([0.5]*mm*nn*2)
        error0 = shape_difference_nondimensional(x0, n, m, original_data)

        def function(x):
            return(shape_difference_nondimensional(x, n, m, original_data) /
                   error0)
        solution = minimize(function, x0, bounds=bounds,
                            options={'maxfun': 30000, 'eps': 1e-07})
        x = solution['x']
        error = solution['fun']

    B = _process_input(x, m=m, n=n, case='coefficients')
    return B, error


def shape_difference_nondimensional(x, n, m, original_data):
    '''Calculates the least square error for TE AND LE.
       - x: [twist,  chord, sweep, shear, x_twist]'''

    # Calculate modeled shape
    [Bu, Bl] = _process_input(x, m=m, n=n, case='coefficients')
    current = CST_3D({'upper': Bu, 'lower': Bl}, mesh=(10, 10))
    current = np.vstack([current['lower'], current['upper']])

    # Calculate error
    [d, i_original, i_model] = directed_hausdorff(original_data, current)
    print(d)
    return d


if __name__ == '__main__':
    # Determining edges
    directory = "./examples/"
    filename = directory + "edges.p"
    # partname = 'wing_right'
    # Read LE and TE
    [Raw_LE, Raw_TE] = pickle.load(open(filename, "rb"), encoding='latin1')
    # Raw_LE = Raw_LE[partname]
    # Raw_TE = Raw_TE[partname]

    # Known values
    cp = ControlPoints()
    cp.half_span = 1.
    cp.eta = [0., 1.0]
    cp.N1 = [0.5, 0.5]
    cp.N2 = [1.0, 1.0]

    # Solve
    [cp, error] = find_ControlPoints(cp, Raw_LE, Raw_TE)

    print('eta', cp.eta)
    print('twist', cp.twist)
    print('chord', cp.chord)
    print('sweep', cp.sweep)
    print('shear', cp.shear)
    print('x_twist', cp.twist_origin['x'])

    pickle.dump(cp, open("edge_properties.p", "wb"))

    [LE_model, TE_model] = edges(cp)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Determining coefficients
    filename = directory + "CST_example.p"
    Data_wing = pickle.load(open(filename, "rb"), encoding='latin1')
    Data_wing = np.vstack([Data_wing['upper'], Data_wing['lower']])

    Data_nondimensional = dimensionalize(Data_wing.copy(), cp, inverse=True)
    # Find non-dimensional domain

    [B, error] = find_coefficients(5, 1, Data_nondimensional)
    print('Bu', B[0])
    print('Bl', B[1])
    Data_CST = CST_3D({'upper': B[0], 'lower': B[1]}, mesh=(10, 10))
    Data_CST = np.vstack([Data_CST['upper'], Data_CST['lower']])

    def np_plot(data, color, plot_type='scatter'):
        x, y, z = data.T
        if plot_type == 'scatter':
            ax1.scatter(x, y, z, c=color)
        elif plot_type == 'plot':
            ax1.plot(x, y, z, c=color)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    # fig2 = plt.figure()
    # ax2 = fig1.add_subplot(111)
    x, y, z = LE_model.T

    np_plot(Raw_LE, 'r')
    np_plot(Raw_TE, 'r')
    np_plot(Data_wing, 'b')
    np_plot(Data_nondimensional, 'k')
    # np_plot(Data_CST, 'm')
    np_plot(LE_model, 'g', 'plot')
    np_plot(TE_model, 'g', 'plot')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()
