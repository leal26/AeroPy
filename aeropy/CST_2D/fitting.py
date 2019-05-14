import math
import numpy as np

from aeropy.geometry.airfoil import CST
from aeropy.xfoil_module import output_reader

from scipy.optimize import fsolve, minimize, differential_evolution


def fitting_shape_coefficients(filename, bounds='Default', n=5,
                               return_data=False, return_error=False,
                               optimize_deltaz=False, solver='gradient',
                               deltaz=None, objective='hausdorf',
                               surface='both', x0=None):
    """Fit shape parameters to given data points
        Inputs:
        - filename: name of the file where the original data is
        - bounds: bounds for the shape parameters. If not defined,
                    Default values are used.
        - n: order of the Bernstein polynomial. If bounds is default
                this input will define the order of the polynomial.
                Otherwise the length of bounds (minus one) is taken into
                consideration"""

    from optimization_tools.hausdorff_distance import hausdorff_distance_2D

    def shape_difference(inputs, optimize_deltaz=False, surface=surface):
        # Define deltaz
        if optimize_deltaz is True or optimize_deltaz == [True]:
            deltasz = inputs[-1]/2.
        else:
            deltasz = deltaz/2.

        # Calculate upper and lower surface
        if surface == 'both':
            y_u = CST(upper['x'], 1, deltasz=deltasz,
                      Au=list(inputs[:n+1]))
            y_l = CST(lower['x'], 1, deltasz=deltasz,
                      Al=list(inputs[n+1:-1]))
        elif surface == 'upper':
            y_u = CST(upper['x'], 1, deltasz=deltasz,
                      Au=list(inputs[:n+1]))
        elif surface == 'lower':
            y_l = CST(lower['x'], 1, deltasz=deltasz,
                      Al=list(inputs[:n+1]))

        # Vector to be compared with
        error = 0
        if surface == 'upper' or surface == 'both':
            a_u = {'x': upper['x'], 'y': y_u}
            if objective == 'hausdorf':
                error += hausdorff_distance_2D(a_u, upper)
            elif objective == 'squared_mean':
                error += np.mean((np.array(a_u['x'])-np.array(upper['x']))**2 +
                                 (np.array(a_u['y'])-np.array(upper['y']))**2)

        if surface == 'lower' or surface == 'both':
            a_l = {'x': lower['x'], 'y': y_l}
            if objective == 'hausdorf':
                error += hausdorff_distance_2D(a_l, lower)
            elif objective == 'squared_mean':
                error += np.mean((np.array(a_l['x'])-np.array(lower['x']))**2 +
                                 (np.array(a_l['y'])-np.array(lower['y']))**2)

        # plt.figure()
        # plt.scatter(a_u['x'], a_u['y'], c='k')
        # plt.scatter(a_l['x'], a_l['y'], c='b')
        # plt.scatter(upper['x'], upper['y'], c='r')
        # plt.scatter(lower['x'], lower['y'], c='g')
        # plt.show()
        return error

    def separate_upper_lower(data):
        for key in data:
            data[key] = np.array(data[key])

        index = np.where(data['y'] > 0)
        upper = {'x': data['x'][index],
                 'y': data['y'][index]}
        index = np.where(data['y'] <= 0)
        lower = {'x': data['x'][index],
                 'y': data['y'][index]}
        # x = data['x']
        # y = data['y']
        # for i in range(len(x)):
        #     if data['y'][i] < 0:
        #         break
        # upper = {'x': x[0:i],
        #          'y': y[0:i]}
        # lower = {'x': x[i:],
        #          'y': y[i:]}
        return upper, lower

    def _determine_bounds_x0(n, optimize_deltaz, bounds):
        if bounds == 'Default':
            upper_bounds = [[0, 1]] + [[-1., 1.]]*n
            lower_bounds = [[0, 1]] + [[-1., 1.]]*n

        if optimize_deltaz:
            if surface == 'both':
                bounds = upper_bounds + lower_bounds + [[0, 0.1]]
                x0 = (n+1)*[0., ] + (n+1)*[0., ] + [0.]
            elif surface == 'upper':
                bounds = upper_bounds + [[0, 0.1]]
                x0 = (n+1)*[0., ] + [0.]
            elif surface == 'lower':
                bounds = lower_bounds + [[0, 0.1]]
                x0 = (n+1)*[0., ] + [0.]
        else:
            bounds = upper_bounds + lower_bounds
            x0 = (n+1)*[0., ] + (n+1)*[0., ]
        return x0, bounds

    # Order of Bernstein polynomial
    if bounds != 'Default':
        n = len(bounds) - 1

    # Obtaining data
    if type(filename) != str:
        data = filename
        x, z = data.T
    elif filename[-2:] == '.p':
        import pickle
        data = pickle.load(open(filename, "rb"), encoding='latin1')
        data = data['wing'][list(data['wing'].keys())[3]]
        x, y, z = data.T
    else:
        data = output_reader(filename, separator='\t', header=['x', 'z'])
        x = data['x']
        z = data['z']

    # Rotating airfoil
    x_TE = (x[0] + x[-1])/2.
    y_TE = (z[0] + z[-1])/2.

    theta_TE = math.atan(-y_TE/x_TE)

    # position trailing edge at the x-axis
    processed_data = {'x': [], 'y': []}
    for i in range(len(x)):
        x_i = x[i]
        z_i = z[i]
        c_theta = math.cos(theta_TE)
        s_theta = math.sin(theta_TE)
        x_rotated = c_theta*x_i - s_theta*z_i
        z_rotated = s_theta*x_i + c_theta*z_i
        processed_data['x'].append(x_rotated)
        processed_data['y'].append(z_rotated)
    data = processed_data

    # determine what is the leading edge and the rotation angle beta
    processed_data = {'x': [], 'y': []}

    min_x = min(x)
    # min_index = data['x'].index(min_x)
    # min_y = data['y'][min_index]

    chord = max(x) - min(x)
    # beta = math.atan((y_TE - min_y)/(x_TE - min_x))

    for i in range(len(x)):
        processed_data['x'].append((x[i] - min_x)/chord)
        processed_data['y'].append(z[i]/chord)
    data = processed_data

    # Determining default bounds
    x0_default, bounds = _determine_bounds_x0(n, optimize_deltaz, bounds)
    if x0 is None:
        x0 = x0_default

    if not optimize_deltaz and deltaz is None:
        deltaz = (data['y'][0] - data['y'][-1])

    if surface == 'both':
        upper, lower = separate_upper_lower(data)
    elif surface == 'upper':
        upper = data
    elif surface == 'lower':
        lower = data

    # Calculate original error
    error0 = shape_difference(x0, optimize_deltaz=optimize_deltaz,
                              surface=surface)

    def f(x):
        return shape_difference(x, optimize_deltaz=optimize_deltaz,
                                surface=surface)/error0
    # Optimize
    if solver == 'differential_evolution':

        result = differential_evolution(f, bounds,
                                        disp=True, popsize=10)
        x = result.x
        f = result.fun
    elif solver == 'gradient':

        solution = minimize(f, x0, bounds=bounds,
                            options={'maxfun': 30000, 'eps': 1e-02})
        x = solution['x']
        f = solution['fun']
    print('order %i  done' % n)

    # Unpackage data
    if surface == 'both' or surface == 'upper':
        Au = list(x[:n+1])
    if surface == 'both':
        if optimize_deltaz:
            Al = list(x[n+1:-1])
            deltaz = x[-1]
        else:
            Al = list(x[n+1:])
    elif surface == 'lower':
        Al = list(x[:n+1])

    # Return Al, Au, and others
    to_return = []
    if return_data:
        to_return.append(data)
    if return_error:
        to_return.append(f)
    to_return.append(deltaz)

    if surface == 'lower' or surface == 'both':
        to_return.append(Al)
    if surface == 'upper' or surface == 'both':
        to_return.append(Au)
    return to_return


def shape_parameter_study(filename, n=5, solver='gradient', deltaz=None,
                          objective='hausdorf', surface='both'):
    """Analyze the shape difference for different Bernstein order
       polynomials.
       - filename: name of dataset to compare with
       - n: Maximum Bernstein polynomial order """
    import pickle

    if deltaz is None:
        optimize_deltaz = True
    else:
        optimize_deltaz = False
    Data = {'error': [], 'Al': [], 'Au': [], 'order': [], 'deltaz': []}
    for i in range(1, n+1):
        if surface == 'both':
            error, deltaz, Al, Au = fitting_shape_coefficients(
                filename, n=i, return_error=True,
                optimize_deltaz=optimize_deltaz, solver=solver, deltaz=deltaz,
                objective=objective, surface=surface)
        elif surface == 'lower':
            error, deltaz, Al = fitting_shape_coefficients(
                filename, n=i, return_error=True,
                optimize_deltaz=optimize_deltaz, solver=solver, deltaz=deltaz,
                objective=objective, surface=surface)
        elif surface == 'upper':
            error, deltaz, Au = fitting_shape_coefficients(
                filename, n=i, return_error=True,
                optimize_deltaz=optimize_deltaz, solver=solver, deltaz=deltaz,
                objective=objective, surface=surface)
        print(error)
        Data['error'].append(error)
        if surface == 'both' or surface == 'lower':
            Data['Al'].append(Al)
        if surface == 'both' or surface == 'upper':
            Data['Au'].append(Au)
        Data['deltaz'].append(deltaz)
        Data['order'].append(i)

    file = open('shape_study.p', 'wb')
    pickle.dump(Data, file)
    return Data
