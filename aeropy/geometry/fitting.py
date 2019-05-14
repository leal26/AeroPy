import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
from scipy.optimize import minimize, differential_evolution
from scipy.spatial.distance import directed_hausdorff
from multiprocessing import Pool
import time

from aeropy.xfoil_module import output_reader


class fitting():
    def __init__(self, **params):
        def _callback(*args):
            return params.get('callback', None)(self.object, *args)
        self.object = params.get('object', np.linspace(10, 50, 9))
        self.update = params.get('update', np.linspace(10, 50, 9))
        self.x0 = params.get('x0', np.linspace(10, 50, 9))
        self.p1 = params.get('p1', np.linspace(10, 50, 9))
        self.p2 = params.get('p2', np.linspace(10, 50, 9))
        self.p1_name = params.get('p1_name', 'Parameter 1')
        self.p2_name = params.get('p2_name', 'Parameter 2')
        self.calculate_points = params.get('calculate_points', 0)
        self.raw = params.get('raw', 0)
        self.callback = _callback

    def convergence_study(self, parallel=True):
        P1_f, P2_f = self._format_parameters()
        if parallel:
            p = Pool(len(P1_f))
        else:
            p = Pool(1)

        input = np.vstack([P1_f, P2_f]).T
        self.solutions = p.map(self.find, input)
        self.error = np.array([self.solutions[i]['fun'] for i in
                               range(len(self.solutions))])
        self.error = self.error.reshape(self.P1.shape)
        self.rel_error = self.error/self.error[0][0]

    def find(self, param_i=[None, None]):
        '''
        inputs: [location, XYZ, sy, ny, xshear]'''
        p1_i, p2_i = param_i
        x0 = self.x0(p1_i, p2_i)
        start = time.time()
        solution = minimize(self.shape_difference, x0, args=param_i)
        end = time.time()
        error = solution['fun']
        if p1_i is None:
            print('error=%f\t time=%f' % (error, end-start))
        else:
            print('p1=%i\t p2=%i\t error=%f\t time=%f' % (p1_i, p2_i, error,
                                                          end-start))
        if self.callback is not None:
            if p1_i is None:
                self.callback()
            else:
                self.callback(p1_i, p2_i)
        return solution

    def shape_difference(self, x, param_i):
        p1_i, p2_i = param_i
        self.update(self.object, x, p1_i, p2_i)
        points = self.calculate_points(self.object, self.raw)

        if self.raw.ndim == 3:
            error = 0
            N = 0
            for i in range(len(points)):
                error += np.linalg.norm(points[i] - self.raw[i])**2
                N += len(points[i])
        else:
            error = np.linalg.norm(points[0] - self.raw)**2
            N = len(points)
        return(error/N)

    def plot_study(self, relative=True):
        if relative:
            z = self.rel_error
        else:
            z = self.error
        fig, ax = plt.subplots()
        cs = ax.contourf(self.P1, self.P2, z, np.linspace(0, 1, 101))

        fig.colorbar(cs, ticks=np.linspace(0, 1, 6))
        # plt.clim(0, 1)
        plt.xlabel(self.p1_name)
        plt.ylabel(self.p2_name)
        plt.show()

    def plot_fit(self):
        network = self.calculate_points(self.object, self.raw)
        plt.figure()
        if network.ndim == 3:
            for i in range(len(network)):
                plt.scatter(network[i, :, 0], network[i, :, 2], c='b',
                            label='fit')
                plt.scatter(self.raw[i, :, 0], self.raw[i, :, 2], c='r',
                            label='raw')
        else:
            plt.scatter(network[:, 0], network[:, 2], c='b', label='fit')
            plt.scatter(self.raw[:, 0], self.raw[:, 2], c='r', label='raw')
        plt.legend()
        plt.show()

    def _format_parameters(self, array_type=int):
        self.p1 = self.p1.astype(array_type)
        self.p2 = self.p2.astype(array_type)
        self.P1, self.P2 = np.meshgrid(self.p1, self.p2)
        P1_f = self.P1.flatten()
        P2_f = self.P2.flatten()
        return(P1_f, P2_f)
