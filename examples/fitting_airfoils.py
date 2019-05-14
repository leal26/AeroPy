from scipy.optimize import minimize
from scipy.spatial.distance import directed_hausdorff
import matplotlib.pyplot as plt
import numpy as np
import math

from aeropy.CST_2D import CST


def fitting_shape_coefficients(data, n=5, upper=True):
    """Fit shape parameters to given data points
        Inputs:
        - filename: name of the file where the original data is
        - bounds: bounds for the shape parameters. If not defined,
                    Default values are used.
        - n: order of the Bernstein polynomial. If bounds is default
                this input will define the order of the polynomial.
                Otherwise the length of bounds (minus one) is taken into
                consideration"""

    def shape_difference(inputs, upper=True):
        if upper:
            y = CST(data[:, 0], 1, deltasz=0, Au=list(inputs))
        else:
            y = CST(data[:, 0], 1, deltasz=0, Al=list(inputs))
        # Vector to be compared with
        a = np.vstack([data[:, 0], np.array(y)]).T

        return np.sqrt(np.mean((a-data)**2))
        # return directed_hausdorff(a, data)[0]

    lower_bounds = [[.001, .2]] + [[-.2, .2]]*(n-1) + [[.001, .2]]

    x0 = .05*np.ones(n+1)

    solution = minimize(shape_difference, x0, bounds=lower_bounds, args=upper)
    Al = solution['x']
    print('order %i  done' % n)

    # Return Al, Au, and others
    return(Al)


def process_data(data):
    # Rotating airfoil
    i_TE = np.argmax(data[:, 0])
    i_LE = np.argmin(data[:, 0])

    data = data - data[i_LE]
    theta = math.atan(-data[i_TE, 1]/data[i_TE, 0])

    temp = data.copy()
    c_theta = math.cos(theta)
    s_theta = math.sin(theta)
    temp[:, 0] = c_theta*data[:, 0] - s_theta*data[:, 1]
    temp[:, 1] = s_theta*data[:, 0] + c_theta*data[:, 1]
    data = temp.copy()

    chord = data[i_TE, 0] - data[i_LE, 0]
    data = (data - data[i_LE])/chord

    data = data[data[:, 0].argsort()]
    return(data, chord, theta)


def shape_parameter_study(filename, n=5):
    """Analyze the shape difference for different Bernstein order
       polynomials.
       - filename: name of dataset to compare with
       - n: Maximum Bernstein polynomial order """
    import pickle

    Data = {'error': [], 'Al': [], 'Au': [], 'order': [], 'deltaz': []}
    for i in range(1, n+1):
        error, deltaz, Al, Au = fitting_shape_coefficients(
            filename, n=i, return_error=True, optimize_deltaz=True)
        Data['error'].append(error)
        Data['Al'].append(Al)
        Data['Au'].append(Au)
        Data['deltaz'].append(deltaz)
        Data['order'].append(i)

    file = open('shape_study.p', 'wb')
    pickle.dump(Data, file)
    return Data


if __name__ == "__main__":
    directory = 'D:\\GitHub\\AeroPy\\examples\\'
    filename = 'naca641212_upper.txt'
    data_upper = np.genfromtxt(directory + filename)
    filename = 'naca641212_lower.txt'
    data_lower = np.genfromtxt(directory + filename)
    Al = fitting_shape_coefficients(data_lower, n=5, upper=False)
    Au = fitting_shape_coefficients(data_upper, n=5, upper=True)
    print(Au)
    print(Al)
    y_u = CST(data_upper[:, 0], 1, deltasz=0, Au=Au)
    y_l = CST(data_lower[:, 0], 1, deltasz=0, Al=Al)
    plt.figure()
    plt.scatter(data_upper[:, 0], data_upper[:, 1], label='raw_upper')
    plt.scatter(data_lower[:, 0], data_lower[:, 1], label='raw_lower')
    plt.plot(data_upper[:, 0], y_u, label='upper')
    plt.plot(data_lower[:, 0], y_l, label='lower')
    plt.legend()
    plt.show()

    # directory = 'D:\\GitHub\\AeroPy\\examples\\JAXA_files\\raw\\'
    # filename = 'airfoil'
    # for i in range(1):
    #     print(i)
    #     data = np.genfromtxt(directory+filename+'_%i.csv' % i, delimiter=',')
    #     data, chord, theta = process_data(data)
    #     Al, Au = fitting_shape_coefficients(data, n=5)
    #     print(Au)
    #     print(Al)
    #     y_u = CST(data[:, 0], 1, deltasz=0, Au=Au)
    #     y_l = CST(data[:, 0], 1, deltasz=0, Al=Al)
    #     plt.figure()
    #     plt.scatter(data[:, 0], data[:, 1], label='raw')
    #     plt.plot(data[:, 0], y_u, label='upper')
    #     plt.plot(data[:, 0], y_l, label='lowerer')
    #     plt.legend()
    #     plt.show()
