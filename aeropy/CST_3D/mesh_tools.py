"""Turns a surface description into a panair network"""
import numpy as np
import sys
from math import sin, cos, sqrt


def axisymmetric_surf(data_x, data_r, N_theta):

    theta_start = np.pi  # radians
    theta_end = np.pi/2.  # radians

    data_t = np.linspace(theta_start, theta_end, N_theta)

    surf_coords = np.zeros([len(data_t), len(data_x), 3])

    for i, t in enumerate(data_t):
        for j, x in enumerate(data_x):
            surf_coords[i, j, 0] = x
            surf_coords[i, j, 1] = data_r[j]*sin(t)
            surf_coords[i, j, 2] = data_r[j]*cos(t)

    num_points = len(data_x)

    max_axial = 200
    num_network = int(num_points/max_axial)
    if not (num_points % max_axial) == 0:
        num_network += 1
    nn = int(num_points/num_network)

    network_list = []
    if num_network > 1:
        for i in range(num_network):
            if i == num_network-1:
                network_list.append(surf_coords[:, i*nn:])
            else:
                network_list.append(surf_coords[:, i*nn:(i+1)*nn+1])

    else:
        network_list.append(surf_coords)

    return network_list


def generate_wake(te_points, x_end, n_points=10, angle_of_attack=0.,
                  user_spacing=None):
    # check that x_end is downstream of all trailing edge points
    if not np.all(te_points[:, 0] < x_end):
        raise RuntimeError("wake must terminate downstream of trailing edge")

    if user_spacing is None:
        spacing = np.linspace
    else:
        spacing = user_spacing

    Ny = te_points.shape[0]
    wake = np.zeros((n_points, Ny, 3))
    aoa_r = angle_of_attack*np.pi/180.
    for j, p in enumerate(te_points):
        x_te, y_te, z_te = p
        length = (x_end-x_te)/np.cos(aoa_r)
        X_0 = spacing(0., length, n_points)
        X_r = X_0*np.cos(aoa_r)
        Z_r = X_0*np.sin(aoa_r)
        wake[:, j, 0] = x_te+X_r
        wake[:, j, 1] = y_te
        wake[:, j, 2] = z_te+Z_r

    return wake


def constant_eta_edge(eta, n_points):  # , cos_spacing=True):
    edge = np.full((n_points, 2), eta)
    cos_space = cosine_spacing()
    edge[:, 0] = cos_space(0., 1., n_points)

    return edge


def constant_psi_edge(psi, n_points, eta_lim):  # , cos_spacing=True):
    edge = np.full((n_points, 2), psi)
    eta0, eta1 = eta_lim
    cos_space = cosine_spacing()
    edge[:, 1] = cos_space(eta0, eta1, n_points)

    return edge


def gen_network_edge(surface, intersection, N_b=None, N_t=None,
                     vertical=False):
    """takes in intersection(s) and returns properly ordered and
    complete network edge

    Parameters
    ----------
    surface : CST3D
        CST surface corresponding to network

    intersection : tuple of 3 1D numpy arrays
        x, y, and z coordinates of points that define intersection

    N_b : int
        Number of points to use in back of intersection if intersection
        doesn't fully define edge.

    N_f : int
        Number of points to use in front of intersection if intersection
        doesn't fully define edge.

    vertical : bool
        If vertical is false, then the edge runs parallel(ish) to the
        psi axis. In other words, eta as a function of psi. If vertical
        is true, than the edge corresponds to the opposite.

    Returns
    -------
    edge : 2D numpy array
        Points that define network edge

    """
    # Bring intersections into parameter space
    x_i, y_i, z_i = intersection
    psi_i, eta_i = surface.inverse(x_i, y_i, z_i)

    if vertical:
        eta_i, psi_i = _process_edge(eta_i, psi_i, N_b, N_t)
    else:
        psi_i, eta_i = _process_edge(psi_i, eta_i, N_b, N_t)

    edge = np.array([psi_i, eta_i]).T

    return edge


def _process_edge(p1, p2, N_b, N_t):
    # sort by p1
    i_sort = np.argsort(p1)
    p1 = p1[i_sort]
    p2 = p2[i_sort]

    # verify completeness. If close, make complete.
    # Completeness means edge completely cuts through domain.
    if abs(p1[0]-0.) < 1e-10:
        p1[0] = 0.
    elif N_b:
        # p1_b = cosine_spacing(0., p1[0], N_b)
        p1_b = np.linspace(0., p1[0], N_b)
        p2_b = np.full(N_b, p2[0])
        p1 = np.concatenate((p1_b[:-1], p1))
        p2 = np.concatenate((p2_b[:-1], p2))
    else:
        raise RuntimeError("edge not complete")
    if abs(p1[-1]-1.) < 1e-10:
        p1[-1] = 1.
    elif N_t:
        # p1_t = cosine_spacing(p1[-1], 1., N_t)
        p1_t = np.linspace(p1[-1], 1., N_t)
        p2_t = np.full(N_t, p2[-1])
        p1 = np.concatenate((p1, p1_t[1:]))
        p2 = np.concatenate((p2, p2_t[1:]))
    else:
        raise RuntimeError("edge not complete")

    return p1, p2


class _uniform_spacing:
    def __init__(self, limits, index):
        self._limits = limits
        self._i = index
        # print(self._limits)
        # print(self._i)

    def __call__(self, *dummy):
        lower, upper = self._limits
        if lower is not None:
            return lower[:, self._i]
        elif upper is not None:
            return upper[:, self._i]
        else:
            raise RuntimeError("must specify edge to use uniform spacing")


def meshparameterspace(shape=(20, 20), psi_limits=(None, None),
                       eta_limits=(None, None),
                       psi_spacing="linear",
                       eta_spacing="linear",
                       user_spacing=(None, None)):
    """Builds curvilinear mesh inside parameter space.

    :param psi_spacing and eta_spacing:
           - 'linear': uniform spacing on interior of the surface
           - 'cosine': cosine spacing
           - 'uniform': spacing matches the spacing along edge
           - 'user': user spacing that is passed in through user_spacing

    :param psi_limits and eta_limits: only define if 'uniform'. Should be
           points where intersection is located.
    """
    if psi_spacing == "cosine":
        x_spacing = cosine_spacing()
    elif psi_spacing == "linear":
        x_spacing = np.linspace
    elif psi_spacing == "uniform":
        x_spacing = _uniform_spacing(eta_limits, 0)
    elif psi_spacing == "user":
        if user_spacing[0] is not None:
            x_spacing = user_spacing[0]
        else:
            raise RuntimeError("must provide user_spacing w/ psi_spacing=user")
    else:
        raise RuntimeError("specified spacing not recognized")

    if eta_spacing == "cosine":
        y_spacing = cosine_spacing()
    elif eta_spacing == "linear":
        y_spacing = np.linspace
    elif eta_spacing == "uniform":
        y_spacing = _uniform_spacing(psi_limits, 1)
    elif eta_spacing == "user":
        if user_spacing[1] is not None:
            y_spacing = user_spacing[1]
        else:
            raise RuntimeError("must provide user_spacing w/ psi_spacing=user")
    else:
        raise RuntimeError("specified spacing not recognized")

    n_psi, n_eta = shape
    psi_lower, psi_upper = psi_limits
    eta_lower, eta_upper = eta_limits

    # if limits aren't specified, set lower to 0 and upper to 1
    if psi_lower is None:
        psi_lower = np.full((n_eta, 2), 0.)
        eta_min = eta_lower[0, 1] if eta_lower is not None else 0.
        eta_max = eta_upper[0, 1] if eta_upper is not None else 1.
        psi_lower[:, 1] = y_spacing(eta_min, eta_max, n_eta)
    if psi_upper is None:
        psi_upper = np.full((n_eta, 2), 1.)
        eta_min = eta_lower[-1, 1] if eta_lower is not None else 0.
        eta_max = eta_upper[-1, 1] if eta_upper is not None else 1.
        psi_upper[:, 1] = y_spacing(eta_min, eta_max, n_eta)
    if eta_lower is None:
        eta_lower = np.full((n_psi, 2), 0.)
        psi_min = psi_lower[0, 0] if psi_lower is not None else 0.
        psi_max = psi_upper[0, 0] if psi_upper is not None else 1.
        eta_lower[:, 0] = x_spacing(psi_min, psi_max, n_psi)
    if eta_upper is None:
        eta_upper = np.full((n_psi, 2), 1.)
        psi_min = psi_lower[-1, 0] if psi_lower is not None else 0.
        psi_max = psi_upper[-1, 0] if psi_upper is not None else 1.
        eta_upper[:, 0] = x_spacing(psi_min, psi_max, n_psi)

    grid = mesh_curvilinear(psi_lower, psi_upper, eta_lower, eta_upper,
                            x_spacing, y_spacing)

    # TODO: the following probably belongs outside the scope of this class
    # if flip:
    #     grid = np.flipud(grid)

    return grid[:, :, 0], grid[:, :, 1]


def mesh_curvilinear(x_lower, x_upper, y_lower, y_upper, x_spacing, y_spacing):
    # verify that corner points match
    xlyl = np.allclose(x_lower[0], y_lower[0], atol=1e-13, rtol=0.)
    xlyu = np.allclose(x_lower[-1], y_upper[0], atol=1e-13, rtol=0.)
    xuyl = np.allclose(x_upper[0], y_lower[-1], atol=1e-13, rtol=0.)
    xuyu = np.allclose(x_upper[-1], y_upper[-1], atol=1e-13, rtol=0.)

    if not (xlyl and xlyu and xuyl and xuyu):
        print(xlyl, xlyu, xuyl, xuyu)
        print(x_lower[0]-y_lower[0])
        print(x_lower[-1]-y_upper[0])
        print(x_upper[0]-y_lower[-1])
        print(x_upper[-1]-y_upper[-1])
        raise RuntimeError("corner points do not match")

    n_x = y_lower.shape[0]
    n_y = x_lower.shape[0]

    grid = np.zeros((n_x, n_y, 2))

    # boundary points are set to match limits exactly
    grid[0, :] = x_lower
    grid[-1, :] = x_upper
    grid[:, 0] = y_lower
    grid[:, -1] = y_upper

    # inner points are evenly spaced between corresponding limits in x and y
    for i in range(1, n_x-1):
        grid[i, 1:-1, 1] = y_spacing(y_lower[i, 1], y_upper[i, 1], n_y)[1:-1]
    for j in range(1, n_y-1):
        grid[1:-1, j, 0] = x_spacing(x_lower[j, 0], x_upper[j, 0], n_x)[1:-1]

    return grid


class cosine_spacing:
    """Parametric function for obtaining a cosine distribution of points"""

    def __init__(self, offset=0, period=1.):
        self._offset = offset
        self._period = period

    def __call__(self, start, stop, num=50):
        # calculates the cosine spacing
        p = self._period
        offset = self._offset
        index = np.linspace(0., 1., num)
        # spacing = .5*(1.-np.cos(p*np.pi*(index-offset)))
        spacing = ((np.cos(np.pi*offset)-np.cos(np.pi*(p*index+offset))) /
                   (np.cos(np.pi*offset)-np.cos(np.pi*(p+offset))))

        points = start+spacing*(stop-start)

        return points

# def cosine_spacing(start, stop, num=50, offset=0, period=1.):
#     # calculates the cosine spacing
#     index = np.linspace(0., 1., num)
#     spacing = .5*(1.-np.cos(period*np.pi*(index-offset)))
#     points = start+spacing*(stop-start)
#     return points


def _distance_point_to_line(P1, P2, PQ):
    x0, y0 = PQ
    x1, y1 = P1
    x2, y2 = P2
    dy = y2-y1
    dx = x2-x1

    return abs(dy*x0-dx*y0+x2*y1-y2*x1)/sqrt(dy*dy+dx*dx)


def _calc_error(point_list):
    # calculates error if all points between endpoints of point_list
    # were removed.
    error = 0.
    front = point_list[0]
    back = point_list[-1]
    for i in range(1, len(point_list)-1):
        error += _distance_point_to_line(front, back, point_list[i])

    return error


def _calc_length(point_list):
    # calculates error if all points between endpoints of point_list
    # were removed.
    x_f, y_f = point_list[0]
    x_b, y_b = point_list[-1]

    length = sqrt((x_b-x_f)**2+(y_b-y_f)**2)

    return length


def coarsen_axi(data_x, data_r, tol, max_length):
    # move x and r data into a list of "points"
    point_list = []
    for i in range(len(data_x)):
        point_list.append(np.array([data_x[i], data_r[i]]))

    # ITERATIVE ALGORITHM
    # Indices for the start and end points of the algorithm
    Pstart = 0
    Pend = len(point_list)-1
    # Indices for 2 pointers that define current range being examined
    P1 = Pstart
    P2 = Pstart+2

    new_point_list = [point_list[Pstart]]

    while P2 <= Pend:
        error = _calc_error(point_list[P1:P2+1])

        if error > tol:
            new_point_list.extend(point_list[P1+1:P2+1])
            P1 = P2
            P2 = P1 + 2
        else:
            while error < tol and P2 <= Pend:
                P2 += 1
                error = _calc_error(point_list[P1:P2+1])
                cell_length = _calc_length(point_list[P1:P2+1])
                # print(cell_length)
                if cell_length > max_length:
                    error += tol*10.
            P2 -= 1
            new_point_list.append(point_list[P2])
            P1 = P2
            P2 = P1 + 2

    if not (new_point_list[-1][0] == point_list[-1][0]):
        new_point_list.append(point_list[-1])

    # print("size of new list", len(new_point_list))
    sys.stdout.flush()
    new_x = np.zeros(len(new_point_list))
    new_r = np.zeros(len(new_point_list))
    for i in range(0, len(new_point_list)):
        new_x[i] = new_point_list[i][0]
        new_r[i] = new_point_list[i][1]

    return new_x, new_r
