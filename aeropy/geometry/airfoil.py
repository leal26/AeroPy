# -*- coding: utf-8 -*-
"""
Created on Fri Oct 09 17:31:59 2015

Current functionalities:
- CST
- flaps
- rotate
@author: Pedro Leal
"""
from __future__ import print_function

import math
import numpy as np


def create_x(c, n=230, distribution='linear'):
    """ Create set of points along the chord befitting for Xfoil. The x
    output is conviniently ordered from TE to LE.

    :param c: float value for the chord.
    :param n: number of points
    :param distribution: linear or polar. linear uses a given delta x to
           find the values of x. Polar uses a delta theta to find values
           of theta that are then used in a circle equation to find the
           x points. Usuallly good for airfoils.

    :rtype: x: numpy.array of values of x

    Because Xfoil it is not efficient or sometimes possible to create a
    uniform distribution of points because the front part of the
    airfoil requires a big amount of points to create a smooth surface
    representing the round leading edge. However there is a maximum of
    points that Xfoil will accept, therefore there is a need to
    overcome this obstacle. This was done by dividing the airfoil in
    to 3 parts.

    - tip: correpondent to 0 to .3% of the chord length, it is the
      most densily populated part of the airfoil to compensate the
      wide variation of slopes

    - middle: correpondent to .3% to 30% of the chord length (
      Shortly above a quarter of the chord length). The second most
      densily populated and the segment with the most amount of
      points. Represents the round section of the airfoil except
      for the tip. Such an amount of points is necessary because it
      is where most of the lift is generated.

    - endbody: correspondent to 30% to 100% of the chord length.
      The less densily populated section. Such unrefined mesh is
      possible because of the simplistic geometry of endbody's
      airfoil, just straight lines.

    >>> print create_x(1.0)
        array([  1.00000000e+00,   9.82051282e-01,   9.64102564e-01,
        ...
        6.66666667e-04,   3.33333333e-04,   0.00000000e+00])

    Created on Thu Feb 27 2014
    @author: Pedro Leal
    """

    if distribution == 'linear':
        max_point = c/4.
        limit = max_point + 0.05*c
        nose_tip = 0.003*c

        # Amount of points for each part (this distribution is empirical)
        N_tip = 10*n/230
        N_middle = 180*n/230
        N_endbody = 40*n/230

        x_endbody = np.linspace(c, limit, N_endbody)
        x_middle = np.linspace(limit, nose_tip, N_middle)
        x_tip = np.linspace(nose_tip, 0, N_tip)

        # Organizing the x lists in a unique list without repeating any
        # numbers
        x2 = np.append(x_middle, np.delete(x_tip, 0))
        x = np.append(x_endbody, np.delete(x2, 0))

    elif distribution == 'polar':
        r = c/2.
        x0 = r
        theta = np.linspace(0, math.pi, n)

        x = x0 + r*np.cos(theta)
    return x

# ===========================================================================
# The following functions are related to creating the airfoil outer mold
# ===========================================================================


def CST(x, c, deltasz=None, Au=None, Al=None, N1=0.5, N2=1., deltasLE=None):
    """
    Based on the paper "Fundamental" Parametric Geometry Representations for
    Aircraft Component Shapes" from Brenda M. Kulfan and John E. Bussoletti.
    The code uses a 1st order Bernstein Polynomial for the "Class Function" /
    "Shape Function" airfoil representation.

    The degree of polynomial is dependant on how many values there are for the
    coefficients. The algorithm is able to use Bernstein Polynomials for any
    order automatically. The algorithm is also able to analyze only the top or
    lower surface if desired. It will recognize by the inputs given. i.e.:
    for CST(x=.2,c=1.,deltasx=.2,Au=.7), there is only one value for Au, so it
    is a Bernstein polynomial of 1st order for the upper surface. By ommiting
    Al the code will only input and return the upper surface.

    Although the code is flexible, the inputs need to be coesive. len(deltasz)
    must be equal to the number of surfaces. Au and Al need to have the same
    length if a full analysis is being realized.

    The inputs are:

        - x:list or numpy. array of points along the chord, from TE and the LE,
          or vice-versa. The code works both ways.

        - c: chord

        - deltasz: list of thicknesses on the TE. In case the upper and lower
          surface are being analyzed, the first element in the list is related
          to the upper surface and the second to the lower surface. There are
          two because the CST method treats the airfoil surfaces as two
          different surfaces (upper and lower)

        - Au: list/float of Au coefficients, which are design parameters. If
          None,the surface is not analyzed. len(Au) equals

        - Al: list/float of Al coefficients, which are design parameters. If
          None,the surface is not analyzed. len(Al) equals

        - deltasLE: leading edge thickness in case it is necessary

    The outputs are:
        - y:
          - for a full analysis: disctionary with keys 'u' and 'l' each with
            a list of the y positions for a surface.
          - for a half analysis: a list with the list of the y postions of the
            the desired surface

    Created on Sun Jan 19 16:36:55 2014

    Updated on Mon May 19 18:13:26 2014

    @author: Pedro Leal
    """

    # Bersntein Polynomial
    def K(r, n):
        K = math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K

    # Shape Function
    def S(r, n, psi):
        S = K(r, n)*(psi**r)*(1.-psi)**(n-r)
        return S

    # Class Function
    def C(N1, N2, psi):
        C = ((psi)**N1)*((1.-psi)**N2)
        return C

    if type(x) == list:
        x = np.array(x)
    # Adimensionalizing
    psi = x/c

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Class Function
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # The Coefficients for an airfoil with a rounded leading edge and a sharp
    # trailing edge are N1=0.5 and N2=1.0.
    C = C(N1, N2, psi)

    # ==========================================================================
    #                   Defining the working surfaces
    # ==========================================================================
    deltaz = {}
    deltaLE = {}
    eta = {}
    y = {}
    Shape = {}

    if Al is not None and Au is not None:
        deltaz['u'] = deltasz[0]
        deltaz['l'] = deltasz[1]
        if deltasLE != None:
            deltaLE['u'] = deltasLE[0]
            deltaLE['l'] = deltasLE[1]
        if len(Au) != len(Al):
            raise Exception("Au and Al need to have the same dimensions")
        elif len(deltasz) != 2:
            raise Exception(
                "If both surfaces are being analyzed, two values for deltasz are needed")

    elif Au is not None and Al is None:
        if type(deltasz) == list:
            if len(deltaz['u']) != 1:
                raise Exception(
                    "If only one surface is being analyzed, one value for deltasz is needed")
            else:
                deltaz['u'] = float(deltasz)
                if deltasLE != None:
                    deltaLE['u'] = float(deltasLE[0])
        else:
            deltaz['u'] = deltasz
            if deltasLE != None:
                deltaLE['u'] = deltasLE
    elif Al is not None and Au is None:
        if type(deltasz) == list:
            if (deltaz['l']) != 1:
                raise Exception(
                    "If only one surface is being analyzed, one value for deltasz is needed")
            else:
                deltaz['l'] = float(deltasz)
                if deltasLE != None:
                    deltaLE['l'] = float(deltasLE)
        else:
            deltaz['l'] = deltasz
            if deltasLE != None:
                deltaLE['l'] = deltasLE
    else:
        raise Exception("Au or Al need to have at least one value")
    # In case leading edge thickness is not used, just make it zero and not use it
    if deltasLE == None:
        deltaLE = {'u': 0, 'l': 0}
    A = {'u': Au, 'l': Al}
    for surface in ['u', 'l']:
        if A[surface] is not None:

            if type(A[surface]) == int or type(A[surface]) == float:
                A[surface] = [A[surface]]
            # the degree of the Bernstein polynomial is given by the number of
            # coefficients
            n = len(A[surface])-1
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Shape Function
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Shape[surface] = 0
            for i in range(len(A[surface])):
                #                print A
                #                print S(i,n,psi)
                if surface == 'l':
                    Shape[surface] -= A[surface][i]*S(i, n, psi)
                else:
                    Shape[surface] += A[surface][i]*S(i, n, psi)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Airfoil Shape (eta=z/c)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Airfoil Shape (eta=z/c)
            if surface == 'l':
                eta[surface] = C*Shape[surface]-psi*deltaz[surface]/c - (1.-psi)*deltaLE[surface]/c
            else:
                eta[surface] = C*Shape[surface]+psi*deltaz[surface]/c + (1.-psi)*deltaLE[surface]/c
            # Giving back the dimensions
            y[surface] = c*eta[surface]
    if Al is not None and Au is not None:
        return y
    elif Au is not None:
        return y['u']
    else:
        return y['l']


def Naca00XX(c, t, x_list, TE_t=False, return_dict='y', for_xfoil=True):
    """
    Generates a simetric NACA airfoil.
    Inputs:
    :param c: chord
    :param t: max thickness normalized by the chord (t/c)
    :param x_list: points between 0 and c (chord lenght). If dictionary with keys
                    u and l, calculates x based on that.
    :param TE_t: trailing edge thickness. If False it is the standard
                 thickness for a NACA airfoil. Otherwise it is the given value.
    :param return_dict: return dictionary y('y') or x and y ('xy')
    Returns dictionary y coordinates with keys 'u'pper and 'l'ower

    The Naca function can be found in: https://en.wikipedia.org/wiki/NACA_airfoil

    Created on Wed Feb 03 12:50:52 2016

    @author: Endryws and Pedro Leal
    """
    y_upper = []
    y_lower = []
    if type(x_list) == list:
        for x in x_list:
            xc = x/c  # Is just for increase speed and facilitate future changes.
            a1 = 5*t*c
            t1 = 0.2969*(math.sqrt(xc))
            t2 = -0.1260*xc
            t3 = -0.3516*(xc**2)
            t4 = 0.2843*(xc**3)
            if TE_t == False:
                t5 = -0.1015*(xc**4)
            else:
                t5 = (TE_t/(10*t) - 0.1036)*(xc**4)
            y = (a1*(t1+t2+t3+t4+t5))
            y_upper.append(y)
            y_lower.append(y*(-1))  # is just for pick the y axis
            # negative numbers
    elif type(x_list) == dict:
        for x in x_list['u']:
            xc = x/c  # Is just for increase speed and facilitate future changes.
            a1 = 5*t*c
            t1 = 0.2969*(math.sqrt(xc))
            t2 = -0.1260*xc
            t3 = -0.3516*(xc**2)
            t4 = 0.2843*(xc**3)
            if TE_t == False:
                t5 = -0.1015*(xc**4)
            else:
                t5 = (TE_t/(10*t) - 0.1036)*(xc**4)
            y = (a1*(t1+t2+t3+t4+t5))
            y_upper.append(y)

        for x in x_list['l']:
            xc = x/c  # Is just for increase speed and facilitate future changes.
            a1 = 5*t*c
            t1 = 0.2969*(math.sqrt(xc))
            t2 = -0.1260*xc
            t3 = -0.3516*(xc**2)
            t4 = 0.2843*(xc**3)
            if TE_t == False:
                t5 = -0.1015*(xc**4)
            else:
                t5 = (TE_t/(10*t) - 0.1036)*(xc**4)
            y = (a1*(t1+t2+t3+t4+t5))
            y_lower.append(y*(-1))  # is just for pick the y axis
            # negative numbers
    if len(x_list) == 1:
        y = {'u': y_upper[0],
             'l': y_lower[0]}
    else:
        if for_xfoil:
            y = {'u': y_upper[:-1],
                 'l': y_lower[-2::-1]}
        else:
            y = {'u': y_upper,
                 'l': y_lower}

    if return_dict == 'y':
        return y

    elif return_dict == 'xy':
        if type(x_list) == list:
            if for_xfoil:
                x = {'u': x_list[:-1],
                     'l': x_list[-2::-1]}
            else:
                x = {'u': x_list,
                     'l': x_list}
            return x, y
        elif type(x_list) == dict:
            return x_list, y


def NACA_four_digit(chord, t, p, m):

    # Mean line coordinates
    x_fwd = np.linspace(0, p*chord)
    x_aft = np.linspace(p*chord, chord)

    y_fwd = (m/p**2)*(2*p*(x/chord)-(x/chord)**2)
    y_aft = (m/(1-p)**2)*((1-2*p)+2*p*(x/chord)-(x/chord)**2)

    # Thickness coordinates
    yt_fwd = Naca00XX(chord, t, x_fwd)
    yt_aft = Naca00XX(chord, t, x_fwd)

    # Angle calculation so that thickness is perpendicular
    dydx_fwd = (2*m/p**2)*(p-x/chord)
    dydx_aft = (2*m/(1-p)**2)*(p-x/chord)

    theta_fwd = np.atan(dydc_fwd)
    theta_aft = np.atan(dydc_aft)

    # Calculate airfoil coordinates
    y_upper = {'x': np.append(x_fwd - yt_fwd*np.sin(theta),
                              x_aft - yt_aft*np.sin(theta)),
               'y': np.append(y_fwd + yt_fwd*np.cos(theta),
                              y_aft + yt_aft*np.cos(theta))}
    y_lower = {'x': np.append(x_fwd + yt_fwd*np.sin(theta),
                              x_aft + yt_aft*np.sin(theta)),
               'y': np.append(y_fwd - yt_fwd*np.cos(theta),
                              y_aft - yt_aft*np.cos(theta))}
    return y_upper, y_lower
# ==============================================================================
# The following function are related to the use of plain flaps
# ==============================================================================


def find_hinge(x_hinge, upper, lower):
    """From the points of upper and lower surface find the y coordinate of
    the hinge at x_hinge

    :param x_hinge: float x-coordinate of the hinge
    :param upper: dictionary with keys x and y, coordiantes of upper surface
    :param lower: dictionary with keys x and y, coordiantes of lower surface
    """
    def find_closest(data, x_hinge, fwd, aft):
        """ Find the points afterwards(aft) and forwards(fwd) that are closest
        to x_hinge. Retuns dictionaries aft and fwd, each with keys x and y."""
        for i in range(len(data['x'])):
            xi = data['x'][i]
            yi = data['y'][i]
            error = abs(xi - x_hinge)

            # If equal just save and break it
            if x_hinge == xi:
                aft['x'] = xi
                aft['y'] = yi
                aft['error'] = error
                fwd['x'] = xi
                fwd['y'] = yi
                fwd['error'] = error
                break

            if xi > x_hinge and error < aft['error']:
                aft['x'] = xi
                aft['y'] = yi
                aft['error'] = error

            if xi < x_hinge and error < fwd['error']:
                fwd['x'] = xi
                fwd['y'] = yi
                fwd['error'] = error
        return fwd, aft

    # Find y hinge
    upper_fwd = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    upper_aft = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    lower_fwd = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    lower_aft = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    hinge = {'x': x_hinge, 'y': 9e9}

    upper_fwd, upper_aft = find_closest(upper, x_hinge, upper_fwd, upper_aft)
    lower_fwd, lower_aft = find_closest(lower, x_hinge, lower_fwd, lower_aft)

    # Interpolate
    hinge['y_upper'] = upper_fwd['y'] + (
        hinge['x'] - upper_fwd['x'])*(
        upper_aft['y'] - upper_fwd['y'])/(
        upper_aft['x'] - upper_fwd['x'])

    hinge['y_lower'] = lower_fwd['y'] + (
        hinge['x'] - lower_fwd['x'])*(
        lower_aft['y'] - lower_fwd['y'])/(
        lower_aft['x'] - lower_fwd['x'])
    hinge['y'] = (hinge['y_upper'] + hinge['y_lower'])/2.

    return hinge


def find_flap(data, hinge, extra_points=None):
    """Create the static airfoil and flap dictionaries containing the outer
    mold coordinates of both. Because it is necessary to have an intersection
    between the static lower surface and the flap lower surface, it is
    sometimes interesting to have an extra poin in the static surface to
    garantee the intersection.

    :param data: dictionary with x and y cooridnates of the whole outer mold

    :param hinge: dictionary with x and y coordinates of hinge. Can be found
                  via find_hinge.

    :param extra_points: include extra points to surface, extending it more
                         than it will actually be (usefull for intersecting
                         lines). Avaialble options: 'lower', 'upper' and None"""
    # Generate empty dictionaries which will store final data.
    flap_data = {'x': [], 'y': []}
    static_data = {'x': [], 'y': []}

    type = None

    for i in range(len(data['x'])):
        xi = data['x'][i]
        yi = data['y'][i]
    # Because of the way that xfoil works, the upper list will always
    # begin from the trailing edge
        if xi > hinge['x']:
            if type == None:
                type = 'upper'
            elif type == 'lower':
                flap_data['x'].append(hinge['x'])
                flap_data['y'].append(hinge['y_' + type])
                static_data['x'].append(hinge['x'])
                static_data['y'].append(hinge['y_' + type])
                # this will make the at the hinge to be included only once
                type = 'upper'
                # If the extra point is unecessary in the lower
                if extra_points == 'lower':
                    static_data['x'].append(data['x'][i+1])
                    static_data['y'].append(data['y'][i+1])
                    static_data['x'].append(data['x'][i+2])
                    static_data['y'].append(data['y'][i+2])
            flap_data['x'].append(xi)
            flap_data['y'].append(yi)

    # Because of the way that xfoil works, the lower list will always
    # begin from the  leading edge
        else:
            if type == None:
                type = 'lower'
            elif type == 'upper':
                # If the extra point is unecessary in the upper
                if extra_points == 'upper' and len(static_data['x']) == 0:
                    static_data['x'].append(data['x'][i-2])
                    static_data['y'].append(data['y'][i-2])
                    static_data['x'].append(data['x'][i-1])
                    static_data['y'].append(data['y'][i-1])
                flap_data['x'].append(hinge['x'])
                flap_data['y'].append(hinge['y_' + type])
                static_data['x'].append(hinge['x'])
                static_data['y'].append(hinge['y_' + type])
                # this will make the at the hinge to be included only once
                type = 'lower'

            static_data['x'].append(xi)
            static_data['y'].append(yi)
    return static_data, flap_data


def find_edges(x, y, both_surfaces=False):
    '''Defining chord as the greatest distance from the trailing edge,
       find the leading edge'''

    x = list(x)
    y = list(y)
    # The Trailing edge will always be the point with greatest x for small angles
    if both_surfaces == True:
        TE_x = (x[0]+x[-1])/2
        TE_y = (y[0]+y[-1])/2
    else:
        TE_x = max(x)
        TE_index = x.index(TE_x)
        TE_y = y[TE_index]
    chord = 0
    LE_index = 0

    for i in range(len(x)):
        distance = math.sqrt((x[i]-TE_x)**2+(y[i]-TE_y)**2)
        if distance > chord:
            LE_index = i
            chord = distance
    theta = math.atan2(TE_y - y[LE_index],
                       TE_x - x[LE_index])
    return ({'x': x[LE_index], 'y': y[LE_index]},
            {'x': TE_x, 'y': TE_y},
            theta, chord)


def rotate(upper, lower={'x': [], 'y': []}, origin={'x': 0, 'y': 0},
           theta=None, unit_theta='rad', move_to_origin=False,
           chord=None, both_surfaces=False):
    """
    :param upper: dictionary with keys x and y, each a list

    :param lower: dictionary with heys x and y, each a list

    :param origin: dictionary with keys x and y, each a float

    :param theta: float representing angle in degrees clock-wise

    :param move_to_origin: if true, will establish origin as (0,0)

    :param chord: will normalize results by this chors
    output: rotated_upper, rotated_lower
    """
    # Calculate angle based on trailing edge if angle not defined
    if theta == None:
        origin, TE, theta, chord = find_edges(upper['x'], upper['y'],
                                              both_surfaces)
        # print('theta', theta)
        # theta = - theta
    output = []

    # For trigonometric relations in numpy, theta must be in radians
    if unit_theta == 'deg':
        theta = theta * np.pi/180.
    # Rotation transformation Matrix

    T = [[np.cos(theta), np.sin(theta)],
         [-np.sin(theta), np.cos(theta)]]

    for coordinates in [upper, lower]:
        rotated_coordinates = {'x': [], 'y': []}
        for i in range(len(coordinates['x'])):
            # The rotation must take place with the center of rotation
            # as the origin
            cx = coordinates['x'][i] - origin['x']
            cy = coordinates['y'][i] - origin['y']
            # Rotate
            rot_x = (T[0][0]*cx + T[0][1]*cy)/chord
            rot_y = (T[1][0]*cx + T[1][1]*cy)/chord
            # Store and add back the values of the origin
            if move_to_origin:
                rotated_coordinates['x'].append(rot_x)
                rotated_coordinates['y'].append(rot_y)
            else:
                rotated_coordinates['x'].append(rot_x + origin['x'])
                rotated_coordinates['y'].append(rot_y + origin['y'])
        output.append(rotated_coordinates)
        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.scatter(coordinates['x'],coordinates['y'])
        # plt.show()

    # In case only one surface is of interest

    if len(lower['x']) == 0:
        return output[0]
    else:
        return output[0], output[1]


def clean(upper_static, upper_flap, lower_static, lower_flap, hinge,
          deflection, N=None, return_flap_i=True,
          unit_deflection='rad'):
    """
    Function to remove intersecting lines and to round transition
    on upper surface.

    :param N: number of points on smooth transition. If not defined, the
             gap distance is taken in to consideration to define the
             number of points inserted. If gap is to small, just created
             a node in the middle

    :param return_flap_i: one of the returns is a list with the indexes
                          at which the upper flap surface ends and that
                          the lower flap surface begins.
    :param deflection: flap angular deflection, clockwise positive
    :param unit_deflection: 'rad'ians or 'deg'rees.
    """
    if unit_deflection == 'deg':
        deflection = deflection * np.pi/180.
    if deflection > 0.:
        # The lower surface has a non-smooth transiton so we only have to
        # clean the intersecting lines
        intersection = intersect_curves(lower_static['x'], lower_static['y'],
                                        lower_flap['x'], lower_flap['y'],
                                        input_type='list')
        try:
            intersection_x = intersection[0][0]
            intersection_y = intersection[1][0]
        except:
            import matplotlib.pyplot as plt
            plt.scatter(lower_static['x'], lower_static['y'], c='b')
            plt.scatter(lower_flap['x'], lower_flap['y'], c='r')
            raise Exception('Lower surfaces are not intersecting')

        closest_x = min(lower_flap['x'], key=lambda x: abs(x-intersection_x))
        closest_i = lower_flap['x'].index(closest_x)
        for i in range(len(lower_flap['x'])):
            # From the plots, every point of the flap before the intersection
            # needs to be elimnated.
            if i <= closest_i:
                if (i == closest_i and closest_x < intersection_x) or i != closest_i:
                    lower_flap['x'][i] = None
                    lower_flap['y'][i] = None

        closest_x = min(lower_static['x'], key=lambda x: abs(x-intersection_x))
        closest_i = lower_static['x'].index(closest_x)
        for i in range(len(lower_static['x'])):
            # From the plots, every point of the flap before the intersection
            # needs to be elimnated.
            if i >= closest_i:
                if (i == closest_i and closest_x > intersection_x) or i != closest_i:
                    lower_static['x'][i] = None
                    lower_static['y'][i] = None

        # Eliminatting the None vectors
        lower_flap['x'] = filter(None, lower_flap['x'])
        lower_flap['y'] = filter(None, lower_flap['y'])
        lower_static['x'] = filter(None, lower_static['x'])
        lower_static['y'] = filter(None, lower_static['y'])

        # add intersection points
        lower_static['x'].append(intersection_x)
        lower_static['y'].append(intersection_y)

        R = hinge['y_upper'] - hinge['y']

        upper_smooth = {}
        upper_smooth['x'] = []
        upper_smooth['y'] = []

        # If N is not N, use an adaptive version to it
        chord = min(lower_flap['x']) - min(lower_static['x'])
        # Minimum step for between nodes at the joint
        N_step = 0.002/chord
        if N == None:
            distance = R*deflection
            if distance <= N_step:
                N = 0  # Space too small, just connect points
            else:
                N = int(math.floor(distance/N_step))

        # Need to create points connecting upper part in circle (part of plain flap
        # exposed during flight). The points created are only for the new surface
        # , the connecting points with the other surface have already been
        # created. If N == 0, substitute connecting nodes for an intermediate
        # node

        if N == 0:
            for key in ['x', 'y']:
                upper_flap[key] = upper_flap[key][:-1]
                upper_static[key] = upper_static[key][1:]
            upper_smooth['x'].append(hinge['x'] + R*np.sin(deflection/2.))
            upper_smooth['y'].append(hinge['y'] + R*np.cos(deflection/2.))
        if N != 0:
            for i in range(N):
                theta = ((N-i)/(N+1.)) * deflection
                upper_smooth['x'].append(hinge['x'] + R*np.sin(theta))
                upper_smooth['y'].append(hinge['y'] + R*np.cos(theta))

        # Assembling all together

        modified_airfoil = {'x': [], 'y': []}

        for key in ['x', 'y']:
            modified_airfoil[key] = upper_flap[key] + upper_smooth[key] + \
                upper_static[key] + lower_static[key] + \
                lower_flap[key]
        if return_flap_i == True:
            i = [len(upper_flap[key]) + len(upper_smooth[key]),
                 len(upper_flap[key]) + len(upper_smooth[key]) +
                 len(upper_static[key]) + len(lower_static[key])]
            return modified_airfoil, i
        else:
            return modified_airfoil
    elif deflection < 0.:
        # The upper surface has a non-smooth transiton so we only have to
        # clean the intersecting lines
        intersection = intersect_curves(upper_static['x'], upper_static['y'],
                                        upper_flap['x'], upper_flap['y'],
                                        input_type='list')
        try:
            intersection_x = intersection[0][0]
            intersection_y = intersection[1][0]
        except:
            import matplotlib.pyplot as plt
            plt.scatter(upper_static['x'], upper_static['y'], c='b')
            plt.scatter(upper_flap['x'], upper_flap['y'], c='r')
            raise Exception('Upper surfaces are not intersecting')

        closest_x = min(upper_flap['x'], key=lambda x: abs(x-intersection_x))
        closest_i = upper_flap['x'].index(closest_x)
        for i in range(len(upper_flap['x'])):
            # From the plots, every point of the flap before the intersection
            # needs to be elimnated.
            if i >= closest_i:
                if (i == closest_i and closest_x < intersection_x) or i != closest_i:
                    upper_flap['x'][i] = None
                    upper_flap['y'][i] = None

        closest_x = min(upper_static['x'], key=lambda x: abs(x-intersection_x))
        closest_i = upper_static['x'].index(closest_x)
        for i in range(len(upper_static['x'])):
            # From the plots, every point of the flap before the intersection
            # needs to be elimnated.
            if i <= closest_i:
                if (i == closest_i and closest_x > intersection_x) or i != closest_i:
                    upper_static['x'][i] = None
                    upper_static['y'][i] = None

        # Eliminatting the None vectors
        upper_flap['x'] = filter(None, upper_flap['x'])
        upper_flap['y'] = filter(None, upper_flap['y'])
        upper_static['x'] = filter(None, upper_static['x'])
        upper_static['y'] = filter(None, upper_static['y'])

        # add intersection points
        upper_static['x'] = [intersection_x] + upper_static['x']
        upper_static['y'] = [intersection_y] + upper_static['y']

        R = hinge['y_upper'] - hinge['y']

        lower_smooth = {}
        lower_smooth['x'] = []
        lower_smooth['y'] = []

        # If N is not N, use an adaptive version to it
        chord = min(lower_flap['x']) - min(lower_static['x'])
        # Minimum step for between nodes at the joint
        N_step = 0.002/chord
        if N == None:
            distance = - R*deflection
            if distance <= N_step:
                N = 0  # Space too small, just connect points
            else:
                N = int(math.floor(distance/N_step))
        # Need to create points connecting lower part in circle (part of plain flap
        # exposed during flight). The points created are only for the new surface
        # , the connecting points with the other surface have already been
        # created. If N == 0, substitute connecting nodes for an intermediate
        # node

        if N == 0:
            for key in ['x', 'y']:
                lower_flap[key] = lower_flap[key][1:]
                lower_static[key] = lower_static[key][:-1]
            lower_smooth['x'].append(hinge['x'] + R*np.sin(deflection/2.))
            lower_smooth['y'].append(hinge['y'] - R*np.cos(deflection/2.))

        if N != 0:
            for i in range(N):
                theta = -((i+1)/(N+1.)) * deflection
                lower_smooth['x'].append(hinge['x'] + R*np.sin(theta))
                lower_smooth['y'].append(hinge['y'] - R*np.cos(theta))

        # Assembling all together

        modified_airfoil = {'x': [], 'y': []}
        for key in ['x', 'y']:
            modified_airfoil[key] = upper_flap[key] + upper_static[key] + \
                lower_static[key] + lower_smooth[key] +\
                lower_flap[key]
        if return_flap_i == True:
            i = [len(upper_flap[key]) - 1,
                 len(upper_flap[key]) + len(upper_static[key]) +
                 len(lower_static[key])]
            return modified_airfoil, i
        else:
            return modified_airfoil
    # If deflection equal to zero, just create flap
    elif deflection == 0.:
        # add intersection points
        upper_static['x'] = upper_static['x'][1:]
        upper_static['y'] = upper_static['y'][1:]
        lower_static['x'] = lower_static['x'][:-1]
        lower_static['y'] = lower_static['y'][:-1]

        modified_airfoil = {'x': [], 'y': []}
        for key in ['x', 'y']:
            modified_airfoil[key] = upper_flap[key] + upper_static[key] + \
                lower_static[key] + lower_flap[key]
        if return_flap_i == True:
            i = [len(upper_flap[key]),
                 len(upper_flap[key]) + len(upper_static[key]) +
                 len(lower_static[key])]
            return modified_airfoil, i
        else:
            return modified_airfoil


def intersect_curves(x1, y1, x2, y2, input_type='list',
                     skip_TE_LE=True):
    """
    Find all intersections betweens curves 1 and 2
    :param x1: x data vector for curve 1
    :param y1: y data vector for curve 1
    :param x2: x data vector for curve 1
    :param y2: y data vector for curve 2

    :returns tuple with all intersections.

    source: http://stackoverflow.com/questions/24549247/
    how-to-compute-which-way-data-points-continue-beyond
    -an-intersection
    """
    # convert inputs to numpy arrays
    if input_type == 'list':
        x1 = np.array(x1)
        x2 = np.array(x2)
        y1 = np.array(y1)
        y2 = np.array(y2)
    # number of points in each curve, number of segments is one less,
    # need at least one segment in each curve
    N1 = x1.shape[0]
    N2 = x2.shape[0]

    # get segment presentation (xi, xi+1; xi+1, xi+2; ..)
    xs1 = np.vstack((x1[:-1], x1[1:]))
    ys1 = np.vstack((y1[:-1], y1[1:]))
    xs2 = np.vstack((x2[:-1], x2[1:]))
    ys2 = np.vstack((y2[:-1], y2[1:]))

    # test if bounding-boxes of segments overlap
    mix1 = np.tile(np.amin(xs1, axis=0), (N2-1, 1))
    max1 = np.tile(np.amax(xs1, axis=0), (N2-1, 1))
    miy1 = np.tile(np.amin(ys1, axis=0), (N2-1, 1))
    may1 = np.tile(np.amax(ys1, axis=0), (N2-1, 1))
    mix2 = np.transpose(np.tile(np.amin(xs2, axis=0), (N1-1, 1)))
    max2 = np.transpose(np.tile(np.amax(xs2, axis=0), (N1-1, 1)))
    miy2 = np.transpose(np.tile(np.amin(ys2, axis=0), (N1-1, 1)))
    may2 = np.transpose(np.tile(np.amax(ys2, axis=0), (N1-1, 1)))
    idx = np.where((mix2 <= max1) & (max2 >= mix1) & (miy2 <= may1) &
                   (may2 >= miy1))  # overlapping segment combinations

    # going through all the possible segments
    x0 = []
    y0 = []
    for (i, j) in zip(idx[0], idx[1]):
        # get segment coordinates
        xa = xs1[:, j]
        ya = ys1[:, j]
        xb = xs2[:, i]
        yb = ys2[:, i]
        # ax=b, prepare matrices a and b
        a = np.array([[xa[1] - xa[0], xb[0] - xb[1]], [ya[1] - ya[0], yb[0] - yb[1]]])
        b = np.array([xb[0] - xa[0], yb[0] - ya[0]])
        r, residuals, rank, s = np.linalg.lstsq(a, b, None)
        # if this is not a
        if rank == 2 and not residuals and r[0] >= 0 and r[0] < 1 and r[1] >= 0 and r[1] < 1:
            if r[0] == 0 and r[1] == 0 and i > 0 and j > 0:
                # super special case of one segment point (not the first) in common, need to differentiate between crossing or contact
                angle_a1 = math.atan2(ya[1] - ya[0], xa[1] - xa[0])
                angle_b1 = math.atan2(yb[1] - yb[0], xb[1] - xb[0])

                # get previous segment
                xa2 = xs1[:, j-1]
                ya2 = ys1[:, j-1]
                xb2 = xs2[:, i-1]
                yb2 = ys2[:, i-1]
                angle_a2 = math.atan2(ya2[0] - ya2[1], xa2[0] - xa2[1])
                angle_b2 = math.atan2(yb2[0] - yb2[1], xb2[0] - xb2[1])

                # determine in which order the 4 angle are
                if angle_a2 < angle_a1:
                    h = angle_a1
                    angle_a1 = angle_a2
                    angle_a2 = h
                if (angle_b1 > angle_a1 and angle_b1 < angle_a2 and (angle_b2 < angle_a1 or angle_b2 > angle_a2)) or\
                        ((angle_b1 < angle_a1 or angle_b1 > angle_a2) and angle_b2 > angle_a1 and angle_b2 < angle_a2):
                    # both in or both out, just a contact point
                    x0.append(xa[0])
                    y0.append(ya[0])
            else:
                x0.append(xa[0] + r[0] * (xa[1] - xa[0]))
                y0.append(ya[0] + r[0] * (ya[1] - ya[0]))
    # filter points close to leading edge and trailing edges
    if skip_TE_LE:
        tol = 1e-4
        try:
            x0, y0 = zip(*[(i, j) for i, j in zip(x0, y0) if i < tol and i > 1-tol])
        except:
            x0 = []
            y0 = []

    return (x0, y0)


def separate_upper_lower(x, y, Cp=None, i_separator=None):
    """Return dictionaries with upper and lower keys with respective
    coordiantes. It is assumed the leading edge is frontmost point at
    alpha=0"""
    # TODO: when using list generated by xfoil, there are two points for
    # the leading edge
    def separate(variable_list, i_separator):
        if type(i_separator) == int:
            variable_dictionary = {'upper': variable_list[0:i_separator+1],
                                   'lower': variable_list[i_separator+1:]}
        elif type(i_separator) == list:
            i_upper = i_separator[0]
            i_lower = i_separator[1]

            variable_dictionary = {'upper': variable_list[0:i_upper],
                                   'lower': variable_list[i_lower:]}
        return variable_dictionary
    # If i is not defined, separate upper and lower surface from the
    # leading edge
    if i_separator == None:
        i_separator = x.index(min(x))

    if Cp == None:
        x = separate(x, i_separator)
        y = separate(y, i_separator)
        return x, y
    else:
        x = separate(x, i_separator)
        y = separate(y, i_separator)
        Cp = separate(Cp, i_separator)
        return x, y, Cp


def find_deflection(x_hinge, upper_cruise, lower_cruise,
                    type='Match trailing edge', alpha=0, Reynolds=0,
                    **kwargs):
    """
    Function to calculate the flap deflection.

    Required inputs:
    :param x_hinge
    :param upper_cruise
    :param lower_cruise
    :param type: can be:
                - 'Match trailing edge': will find the deflection where the
                baseline's trailing edge matches that of the objective, i.e.
                match the the trailing edge of a traditional flap with a
                morphing conitunous flap.
                - 'Same Cl': will find deflection where Cl of flapped airfoil
                is equal to the the objective airfoil.
                - 'Range of deflection'
    kwrgs arguments for 'Match trailing edge':
    :param x_TE_objective
    :param y_TE_objective
    :param x_TE_baseline
    :param y_TE_baseline

    kwrgs arguments for 'Same Cl':
    :param Cl_objective

    kwrgs arguments for 'Same Cd':
    :param Cd_objective

    kwrgs arguments for 'Range of deflection'
    :param max_deflection

    Optional arguments:
    :param alpha
    :param Reynolds
    """
    import xfoil_module as xf

    airfoil = 'flapped_airfoil'

    if type == 'Match trailing edge':
        x_TE_objective = kwargs['x_TE_objective']
        y_TE_objective = kwargs['y_TE_objective']
        x_TE_baseline = kwargs['x_TE_baseline']
        y_TE_baseline = kwargs['y_TE_baseline']
    elif type == 'Same Cl':
        Cl_objective = kwargs['Cl_objective']
    elif type == 'Same Cd':
        Cd_objective = kwargs['Cd_objective']
    elif type == 'Range of deflection':
        max_deflection = kwargs['max_deflection']
        init_deflection = kwargs['init_deflection']
        step = kwargs['step']
    hinge = find_hinge(x_hinge, upper_cruise, lower_cruise)

    upper_static, upper_flap = find_flap(upper_cruise, hinge)
    lower_static, lower_flap = find_flap(lower_cruise, hinge)

    # ==============================================================================
    #  Find Deflection
    # ==============================================================================
    if type == 'Match trailing edge':
        # Calculate deflection angle in radians
        deflection = np.arctan2(hinge['y'] - y_TE_objective, x_TE_objective - hinge['x']) - \
            np.arctan2(hinge['y'] - y_TE_baseline, x_TE_baseline - hinge['x'])
        # Convet to degrees
        deflection = deflection*180./np.pi

        upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)

        flapped_airfoil = clean(upper_static, upper_rotated, lower_static,
                                lower_rotated, hinge, deflection, N=5)

        xf.create_input(airfoil, 'Type 3', airfoil=flapped_airfoil)

        xf.call(airfoil, alfas=alpha, Reynolds=Reynolds, iteration=100, GDES=True,
                output='Polar')
        filename = xf.file_name(airfoil, alpha, output='Polar')
        Data = xf.output_reader(filename, output='Polar')

        # Rename to make more convenient
        Cd = Data['CD']
        Cl = Data['CL']

    elif type == 'Same Cl':
        current_Cl = 0
        current_Cd = 0
        previous_Cl = 0
        deflection = 2.
        history = []
        # Rough search
        step = 2.
        while current_Cl < Cl_objective and deflection < 90.:  # Convet to degrees
            previous_Cl = current_Cl
            previous_Cd = current_Cd
            previous_deflection = deflection
            history.append([previous_deflection, previous_Cl, previous_Cd])
            deflection = deflection + step
            print(deflection)
            upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)

            flapped_airfoil = clean(upper_static, upper_rotated, lower_static,
                                    lower_rotated, hinge, deflection, N=5)

            xf.create_input(airfoil, 'Type 3', airfoil=flapped_airfoil)

            xf.call(airfoil, alfas=alpha, Reynolds=Reynolds, iteration=100, GDES=True,
                    output='Polar')
            filename = xf.file_name(airfoil, alpha, output='Polar')
            Data = xf.output_reader(filename, output='Polar')

            # Rename to make more convenient
            current_Cd = Data['CD'][0]
            current_Cl = Data['CL'][0]

        for i in range(1, 6):
            # Fine search
            step = 0.5**i
            current_Cl = previous_Cl
            while current_Cl < Cl_objective and deflection < 90.:  # Convet to degrees
                previous_Cl = current_Cl
                previous_Cd = current_Cd
                previous_deflection = deflection
                history.append([previous_deflection, previous_Cl, previous_Cd])
                deflection = deflection + step

                upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)
                previous_flapped_airfoil = flapped_airfoil
                flapped_airfoil = clean(upper_static, upper_rotated, lower_static,
                                        lower_rotated, hinge, deflection, N=5)

                xf.create_input(airfoil, 'Type 3', airfoil=flapped_airfoil)

                xf.call(airfoil, alfas=alpha, Reynolds=Reynolds, iteration=100, GDES=True,
                        output='Polar')
                filename = xf.file_name(airfoil, alpha, output='Polar')
                Data = xf.output_reader(filename, output='Polar')

                # Rename to make more convenient
                current_Cd = Data['CD'][0]
                current_Cl = Data['CL'][0]
            print(current_Cl, deflection)

    elif type == 'Same Cd':
        current_Cl = 0
        current_Cd = 0
        previous_Cl = 0
        deflection = 2.
        history = []
        # Rough search
        step = 2.
        while current_Cd < Cd_objective and deflection < 90.:  # Convet to degrees
            previous_Cl = current_Cl
            previous_Cd = current_Cd
            previous_deflection = deflection
            history.append([previous_deflection, previous_Cl, previous_Cd])
            deflection = deflection + step
            print(deflection)
            upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)

            flapped_airfoil = clean(upper_static, upper_rotated, lower_static,
                                    lower_rotated, hinge, deflection, N=5)

            xf.create_input(airfoil, 'Type 3', airfoil=flapped_airfoil)

            xf.call(airfoil, alfas=alpha, Reynolds=Reynolds, iteration=100, GDES=True,
                    output='Polar')
            filename = xf.file_name(airfoil, alpha, output='Polar')
            Data = xf.output_reader(filename, output='Polar')

            # Rename to make more convenient
            current_Cd = Data['CD'][0]
            current_Cl = Data['CL'][0]

        for i in range(1, 6):
            # Fine search
            step = 0.5**i
            current_Cd = previous_Cd
            while current_Cd < Cd_objective or deflection < 90.:  # Convet to degrees
                previous_Cl = current_Cl
                previous_Cd = current_Cd
                previous_deflection = deflection
                history.append([previous_deflection, previous_Cl, previous_Cd])
                deflection = deflection + step

                upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)
                previous_flapped_airfoil = flapped_airfoil
                flapped_airfoil = clean(upper_static, upper_rotated, lower_static,
                                        lower_rotated, hinge, deflection, N=5)

                xf.create_input(airfoil, 'Type 3', airfoil=flapped_airfoil)

                xf.call(airfoil, alfas=alpha, Reynolds=Reynolds, iteration=1000, GDES=True,
                        output='Polar')
                filename = xf.file_name(airfoil, alpha, output='Polar')
                Data = xf.output_reader(filename, output='Polar')

                # Rename to make more convenient
                print(deflection)
                current_Cd = Data['CD'][0]
                current_Cl = Data['CL'][0]
            print(current_Cd, deflection)
            print(history)
            print(Cd_objective)
        deflection = previous_deflection
        Cd = previous_Cd
        Cl = previous_Cl
        flapped_airfoil = previous_flapped_airfoil

    elif type == 'Range of deflection':
        current_Cl = 0
        current_Cd = 0
        previous_Cl = 0
        deflection = init_deflection
        history = {'CD': [], 'CL': [], 'deflection': []}
        # Rough search
        while deflection < max_deflection:  # Convet to degrees
            previous_Cl = current_Cl
            previous_Cd = current_Cd
            previous_deflection = deflection

            deflection = deflection + step
            print(deflection)
            upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)

            flapped_airfoil = clean(upper_static, upper_rotated, lower_static,
                                    lower_rotated, hinge, deflection, N=5)

            xf.create_input(airfoil, 'Type 3', airfoil=flapped_airfoil)

            xf.call(airfoil, alfas=alpha, Reynolds=Reynolds, iteration=300, GDES=True,
                    output='Polar')
            filename = xf.file_name(airfoil, alpha, output='Polar')
            Data = xf.output_reader(filename, output='Polar')

            # Rename to make more convenient
            print(deflection)
            current_Cd = Data['CD'][0]
            current_Cl = Data['CL'][0]
            if deflection <= max_deflection:
                history['CD'].append(current_Cd)
                history['CL'].append(current_Cl)
                history['deflection'].append(deflection)
        return history
    return deflection, Cd, Cl, flapped_airfoil


def offset_point(x, y, rho, output_format='separate'):
    """Function to calculate offset curve for a line given by points x and y
    with a distance rho. If output_format = 'separate', the output are
    two list x and y, if 'together, it outputs a unique set with where
    each component is  apari of x and y.'"""
    if type(x) == np.array and type(y) == np.array:
        x_offset = x + rho*y/np.sqrt(x*x + y*y)
        y_offset = y - rho*x/np.sqrt(x*x + y*y)

    if type(x) == list and type(y) == list:
        x_offset = []
        y_offset = []
        for i in range(len(x)):
            # Ignore 0 points because they result in zero denominator
            try:
                x_offset.append(x[i] + rho*y[i]/np.sqrt(x[i]**2 + y[i]**2))
                y_offset.append(y[i] - rho*x[i]/np.sqrt(x[i]**2 + y[i]**2))
            except ZeroDivisionError:
                continue
    # TODO: include option of eliminating lines that intersect

    if output_format == 'together':
        return zip(x_offset, y_offset)
    if output_format == 'separate':
        return x_offset, y_offset

