import math
import numpy as np

from aeropy.geometry.airfoil import CST
from aeropy.xfoil_module import output_reader

from scipy.optimize import fsolve, minimize, differential_evolution
try:
    from optimization_tools.hausdorff_distance import hausdorff_distance_2D
except:
    pass

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
            upper_bounds = [[0, 3]] + [[-3., 3.]]*n
            lower_bounds = [[0, 3]] + [[-3., 3.]]*n

        if optimize_deltaz:
            if surface == 'both':
                bounds = upper_bounds + lower_bounds + [[0, 0.4]]
                x0 = (n+1)*[0., ] + (n+1)*[0., ] + [0.]
            elif surface == 'upper':
                bounds = upper_bounds + [[0, 0.4]]
                x0 = (n+1)*[0., ] + [0.]
            elif surface == 'lower':
                bounds = lower_bounds + [[0, 0.4]]
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
                                        disp=False, popsize=10)
        x = result.x
        f = result.fun
    elif solver == 'gradient':

        solution = minimize(f, x0, bounds=bounds,
                            options={'maxfun': 30000, 'eps': 1e-02})
        x = solution['x']
        f = solution['fun']
    #print('order %i  done' % n)

    # Unpackage data
    if optimize_deltaz:
        deltaz = x[-1]/2.
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


# This is the depricated version of the fitting code
# def fitting_shape_coefficients(data, order, translate=True, rotate=True,
                               # mirror=True, filter_size=10, deltaz=0,
                               # N1=None, N2=None,
                               # solver='differential_evolution',
                               # error='eucledian'):
    # '''Fit shape coefficient to 2D geometry. Assumes data is sorted from
      # TE to LE along upper surface and back to TE along lower surface.
      # Can calculate leading edge for a simple case where upper surface
      # is all positive.

    # :param data: dictionary with x and y keys.

    # :param order: order of the CST equations.
    # '''
    # from scipy.optimize import differential_evolution
    # from optimization_tools import hausdorff_distance_2D, eucledian_shape_difference
    # from scipy.optimize import minimize

    # def total_shape_difference(inputs, error_type='eucledian'):
        # # TODO: update this so that there is an argument upper and lower
        # # from raw data
        # Au = inputs[:order+1]
        # Al = inputs[order+1:2*order+2]

        # if find_class:
            # N1 = inputs[-2]
            # N2 = inputs[-1]

        # else:
            # N1 = 0.5
            # N2 = 1.0

        # error = 0
        # error += shape_difference_upper(inputs, upper, error_type)
        # error += shape_difference_lower(inputs, lower, error_type)

        # return error

    # def shape_difference_upper(inputs, raw_upper, error):
        # shape_coefficients = inputs[:order+1]
        # if find_class:
            # y = CST(raw_upper['x'], 1, deltasz=deltaz/2., Au=shape_coefficients,
                    # N1=inputs[-2], N2=inputs[-1])
        # else:
            # y = CST(raw_upper['x'], 1, deltasz=deltaz/2., Au=shape_coefficients)
        # a = raw_upper
        # b = {'x': raw_upper['x'], 'y': y}

        # if error == 'eucledian':
            # d = eucledian_shape_difference(a, b)
        # else:
            # d = hausdorff_distance_2D(a, b)

        # return d

    # def shape_difference_lower(inputs, raw_lower, error):
        # shape_coefficients = inputs[order+1:2*order+2]

        # if find_class:
            # y = CST(raw_lower['x'], 1, deltasz=deltaz/2., Al=shape_coefficients,
                    # N1=inputs[-2], N2=inputs[-1])
        # else:
            # y = CST(raw_lower['x'], 1, deltasz=deltaz/2., Al=shape_coefficients,)

        # a = raw_lower
        # b = {'x': raw_lower['x'], 'y': y}
        # if error == 'eucledian':
            # d = eucledian_shape_difference(a, b)
        # else:
            # d = hausdorff_distance_2D(a, b)
        # return d

    # def separate_upper_lower(data):
        # for i in range(len(data['x'])):
            # if data['y'][i] < 0:
                # break
        # upper = {'x': data['x'][0:i],
                 # 'y': data['y'][0:i]}
        # lower = {'x': data['x'][i:],
                 # 'y': data['y'][i:]}
        # return upper, lower

    # def processing_data():

        # processed_data = {'x': [], 'y': []}
        # min_x = min(data['x'])
        # chord = max(data['x']) - min(data['x'])

        # # translate to origin and normalize everything
        # for i in range(len(data['x'])):
            # processed_data['x'].append((data['x'][i] - min_x)/chord)
            # processed_data['y'].append(data['y'][i]/chord)
        # data = processed_data

        # # Mirroring data
        # if mirror:
            # processed_data = {'x': [], 'y': []}
            # for i in range(len(data['x'])):
                # processed_data['x'].append(.5 - (data['x'][i] - .5))
                # processed_data['y'].append(data['y'][i])
            # data = processed_data

        # # plt.scatter(data['x'], data['y'])
        # # plt.show()

        # # Tranlating in y
        # if translate:
            # # LE_x = min(data['x'])
            # # LE_index = data['x'].index(LE_x)
            # # LE_y = data['y'][LE_index]
            # for i in range(len(data['x'])):
                # data['y'][i] -= LE_y

        # # plt.scatter(data['x'], data['y'])
        # # plt.show()

        # if rotate:
            # x_TE = (data['x'][0] + data['x'][-1])/2.
            # y_TE = (data['y'][0] + data['y'][-1])/2.

            # # Rotating airfoil
            # theta_TE = math.atan(-y_TE/x_TE)

            # # position trailing edge at the x-axis
            # processed_data = {'x': [], 'y': []}
            # for i in range(len(data['x'])):
                # x = data['x'][i]
                # y = data['y'][i]
                # c_theta = math.cos(theta_TE)
                # s_theta = math.sin(theta_TE)
                # x_rotated = c_theta*x - s_theta*y
                # y_rotated = s_theta*x + c_theta*y
                # processed_data['x'].append(x_rotated)
                # processed_data['y'].append(y_rotated)
            # data = processed_data

        # # Reduce number of points
        # if filter_size != 1:
            # processed_data = {'x': [], 'y': []}
            # for i in range(len(data['x'])):
                # for key in data.keys():
                    # if i % filter_size == 0:
                        # processed_data[key].append(data[key][i])
                    # elif i == len(data['x']) - 1:
                        # processed_data[key].append(data[key][i])
            # data = processed_data

        # for key in data.keys():
            # data[key] = list(reversed(data[key]))

        # deltaz = (data['y'][0] - data['y'][-1])/2.

        # upper, lower = separate_upper_lower(data)
        # return upper, lower

    def find_values(order, solver):
        bounds = [[0., 1.]] + order*[[0., 1.], ] + [[0., 1.]] + order*[[-1., 1.], ]

        # If N1 and N2 are not defined the default values for a subsonic airfoil are
        # used. N1 and N2 are the same for upper and lower.
        if find_class:
            bounds += [[0., 1.], [0., 1.]]
        if solver == 'differential_evolution':
            result = differential_evolution(total_shape_difference, bounds,
                                            disp=True, popsize=40)
            Au = result.x[:order+1]
            Al = result.x[order+1:2*order+2]
            if find_class:
                N1 = result.x[-2]
                N2 = result.x[-1]
            else:
                N1 = 0.5
                N2 = 1.0
            error = solution.fun

        if solver == 'gradient':
            x0 = (2*order+2)*[0.05, ]
            if find_class:
                x0 += [0.5, 1.0]
            solution = minimize(total_shape_difference, x0, bounds=bounds)
            # print solution

            Au = solution['x'][:order+1]
            Al = solution['x'][order+1:2*order+2]
            if find_class:
                N1 = solution['x'][-2]
                N2 = solution['x'][-1]
            else:
                N1 = 0.5
                N2 = 1.0
            error = solution['fun']
        return Au, Al, N1, N2, deltaz, error
    # If data is a list with two dictionaries, it is assumed that the data
    # has been already processed
    if type(data) == list:
        upper = data[0]
        lower = data[1]
    else:
        [upper, lower] = processing_data()

    if N1 is None and N2 is None:
        find_class = True
    # print total_shape_difference([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    # BREAK
    # ==============================================================================
    # Optimizing shape
    # ==============================================================================
    if type(order) == float or type(order) == int:
        [Au, Al, N1, N2, deltaz, error] = find_values(order, solver)
    else:
        Au = []
        Al = []
        N1 = []
        N2 = []

        error = []
        order_list = order
        for i in range(len(order)):
            order = order_list[i]
            [Au_i, Al_i, N1_i, N2_i, deltaz_i, error_i] = find_values(order, solver)
            Au.append(Au_i)
            Al.append(Al_i)
            N1.append(N1_i)
            N2.append(N2_i)
            error.append(error_i)
    return Au, Al, N1, N2, deltaz, error