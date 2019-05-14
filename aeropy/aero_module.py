# -*- coding: utf-8 -*-
"""
Current funcionatilities:
- Lifting line theory
- generate field pressures for Abaqus or other softwares
- air properties calculator
- Reynolds calculator
Created on Mon Jul 20 17:26:19 2015

@author: Pedro Leal
"""
from __future__ import print_function
from __future__ import absolute_import
import math
import numpy as np

#class Wing():
#    def __init__(self, alpha_L_0_root, c_D_xfoil, N=10, b=10., taper=1.,
#                   chord_root=1, alpha_root=0., V=1.):
#        self.alpha_L_0_root = alpha_L_0_root
#        self.c_D_xfoil = c_D_xfoil
#        self.N = N
#        self.b = b
#        self.taper = taper
#        self.chord_root = chord_root
#        self.alpha_root = alpha_root
#        self.V = V

#==============================================================================
# Functions that calculate aerodynamic properties based on already calcualted
# aerodynamic properties from other modules
#==============================================================================
def LLT_calculator(alpha_L_0_root, c_D_xfoil, N=10, b=10., taper=1.,
                   chord_root=1, alpha_root=0., V=1.):
    """
    Calculate the coefficients for a Wing.
    TODO :  - Include elliptical wing
            - When alpha_L_0_root = zero, nan!
            - Include non rectangular wings
            - something else?
            """
    def x_theta_converter(input, b, Output='x'):
        """Converts cartesian coordinate in a polar coordinate."""
        if Output == 'x':
            output = -(b/2) * np.cos(input)
        elif Output == 'theta':
            raise Exception('I did not program this')
        return output

    def geometric_calculator(taper, chord_root, b, N, Type='Linear'):
        """Calculate the following geometric properties:
        - S: Platform Area (currently rectangle)
        - AR: Aspect Ratio
        - c: array containg chord values along all the wing
        """
        theta = np.linspace(np.pi/2, np.pi * N/(N+1), N)
        x = x_theta_converter(theta, b, Output='x')

        c = chord_root * (np.ones(N) - (1-taper)*abs(x)/chord_root)

        if Type == 'Linear':
            S = (1+taper) * (b/2) * chord_root
            AR = b**2/S # Aspect Ratio
        else:
            raise Exception('I did not program this')
        return c, AR

    def A_calculator(N, b, taper, alpha_root, chord_root, alpha_L_0_root):
        """Solve the system of N linear equations with N variables.
        DEFINITION OF A IS IN IPYTHON.
        """
        # Calculate geometric properties
        c, AR = geometric_calculator(taper, chord_root, b, N, Type='Linear')

        # Converting angles to radians
        alpha_root = alpha_root*np.pi/180.
        alpha_L_0_root = alpha_L_0_root*np.pi/180.

        # Avoid using theta = 0,pi,etc, because of zero division
        # Since sine is an odd function and the wing is symmetric, avoid
        # points on the other side of the wing.
        # Otherwise, certain terms would cancel eash other.
        # That is why: N=1+2*j
        theta = np.linspace(np.pi/2, np.pi*N/(N+1), N)
        alpha = alpha_root * np.ones(N)
        alpha_L_0 = alpha_L_0_root * np.ones(N) # CONSTANT AIRFOIL SECTION

        D = np.zeros(N)
        C = np.zeros((N,N))
        for i in range(0,N):
            D[i] = alpha[i] - alpha_L_0[i]
            for j in range(0,N):
                n = 1+2*j
                C[i][j] = ((2*b) / (np.pi*c[i]) + n/np.sin(theta[i])) * \
                          np.sin(n*theta[i])
        A = np.linalg.solve(C, D)

        return A, theta

    def gamma_calculator(b, V, A, theta):
        """Calculate the source strengths."""
        N = len(theta)
        # Calculating gammas
        gamma = []
        for th in theta:
            gamma_temp = 0
            for i in range(0,N):
                n = i+1
                gamma_temp += 2*b*V*A[i]*np.sin(n*th)
            gamma.append(gamma_temp)
        # For tensor manipulation, the data needs to be an np.array object
        gamma = np.array(gamma)
        return gamma

    def coefficient_calculator(A, gamma, c, AR, b, V):
        """Calculate 3D Lift, Drag and efficiency coefficients. The
        section lift coefficient and the lift distribution (roughly equal
        to the pressure distribution) are also caulculated.

        Output: Dictionary with the following keys:
        - cls : section lift coefficient.
        - C_L : 3D Lift Coefficient.
        - C_D : 3D Drag Coefficient.
        - e : efficiency. (as defined in Anderson's Aerodynamics book)
              Between 0 an 1, where 1 is equal to the efficiency of an
              elliptical wing.
        - distribution : distribution along axis perpendicular to the
                         the cross section."""
        N = len(gamma)
        # Calculating section lift coefficients
        cl = 2.*gamma/(c*V)
        # Lift Coefficient
        C_L = A[0]*np.pi*AR
        # Drag Coefficient
        try:
            delta = 0
            for i in range(1,N):
                n = i+1
                delta += n * (A[i]/A[0])**2
            e = 1 / (1+delta)
            C_Di = C_L**2 / (np.pi*e*AR)
        except:
            e = 'Does not make sense. All A coefficients are zero and' \
              ' there is a division by zero at the calculation of delta'
            C_Di = 0

        distribution = cl/cl[0]
        return {'cls':cl, 'C_L':C_L, 'C_Di':C_Di, 'e':e,
                'distribution':distribution}

    def total_drag_calculator(coefficients, c_D_xfoil, b, x):
        """
        From xfoil we have the friction and pressure drag components for
        2D. From LLT we have the 3D induced drag component and the
        pressure distribution. Integrating the 2D over the distribution
        and adding to the 3D induced drag, we obtain the overall drag
        coefficient.

        CURRENTLY WRITTEN FOR SYMMETRIC AIRFOIL."""
        def trapezoid(y,x):
            s = 0
            n = len(y)
            for i in range(0,n-1):
                s += (y[i+1]+y[i]) * (x[i+1]-x[i])/2.
            return s
        # Add the wingtip where the circulation is euqal to zero
        distribution = list(coefficients['distribution'])
        distribution.append(0)
        x = list(x)
        x.append(b/2.)

        C_D_xfoil = 2 * c_D_xfoil * trapezoid(distribution,x)/b
        C_D = C_D_xfoil + coefficients['C_Di']
        return C_D

    # Calculate the Fourrier Constants and theta
    A, theta = A_calculator(N, b, taper, alpha_root, chord_root,
                            alpha_L_0_root)

    # For plotting we calculate the postion along the span
    x = x_theta_converter(theta, b, Output='x')
    # Calculate the circulation coefficients, the gammas
    gamma = gamma_calculator(b, V, A, theta)
    # Calculate certain geometric properties such as chord and AR
    c, AR = geometric_calculator(taper, chord_root, b, N)
    # Calculate the coefficients
    coefficients = coefficient_calculator(A, gamma, c, AR, b, V)
    # Calculate the total drag
    C_D = total_drag_calculator(coefficients, c_D_xfoil, b, x)
    coefficients['C_D'] = C_D
    return coefficients

def calculate_moment_coefficient(x, y, Cp, alpha, c = 1., x_ref = 0.25,
                                 y_ref = 0., flap = False):
    """
    Calculate the moment coeffcient. Inputs are x and y coordinates, and
    pressure coefficients (Cp). Inputs can be in a list in xfoil format
    (counterclockwise starting from the trailing edge, in case necessary,
    check create_input function from xfoil_module) or dictionaries with
    'upper' and 'lower' keys.

    :param flap: if true, also calculates the moment contribution from
                 the trailing edge and the panels in front of the flap
                 (that are not directly in contact with the air)
    """
    def separate_upper_lower(x,y,Cp):
        """Return dictionaries with upper and lower keys with respective
        coordiantes. It is assumed the leading edge is frontmost point at
        alpha=0"""
        #TODO: when using list generated by xfoil, there are two points for
        #the leading edge
        def separate(variable_list, i_separator):
            variable_dictionary = {'upper': variable_list[0:i_separator+1],
                                   'lower': variable_list[i_separator+1:]}
            return variable_dictionary
        i_separator = x.index(min(x))

        x = separate(x, i_separator)
        y = separate(y, i_separator)
        Cp = separate(Cp, i_separator)
        return x, y, Cp
    # If list need to separate in to upper and lower inside a dicitonary
    if type(x) == list and type(y) == list and type(Cp) == list:
        x, y, Cp = separate_upper_lower(x, y, Cp)
    elif type(x) != dict and type(y) != dict and type(Cp) != dict:
        raise Exception("Not all inputs are the same required format (list/dict)")

    #rotating coordinates
    alpha = math.radians(alpha)
    bar_x = {'upper':[], 'lower':[]}
    bar_y = {'upper':[], 'lower':[]}
    for key in x:
        for i in range(len(x[key])):
            bar_x[key].append(x[key][i]*math.cos(alpha) + y[key][i]*math.sin(alpha))
            bar_y[key].append(y[key][i]*math.cos(alpha) - x[key][i]*math.sin(alpha))

    bar_x_ref = x_ref*math.cos(alpha) + y_ref*math.sin(alpha)
    bar_y_ref = y_ref*math.cos(alpha) - x_ref*math.sin(alpha)
    #Rewriting coordinate variables
    x = bar_x
    y = bar_y
    x_ref = bar_x_ref
    y_ref = bar_y_ref


    Cm = 0.
    for key in ['upper', 'lower']:
        for i in range(len(x[key])-1):
            Cm += (1./2*c**2)*(Cp[key][i] + Cp[key][i+1])* \
                  (((x[key][i] + x[key][i+1])/2. - x_ref) *(x[key][i] - x[key][i+1]) + \
                  ((y[key][i] + y[key][i+1])/2. - y_ref) *(y[key][i] - y[key][i+1]))

    if flap == True:
        # Trailing edge contribution
        Cm += (1./2*c**2)*(Cp['upper'][0] + Cp['lower'][-1])* \
              (((x['upper'][0] + x['lower'][-1])/2. - x_ref) *(x['lower'][-1] - x['upper'][0]) + \
               ((y['upper'][0] + y['lower'][-1])/2. - y_ref) *(y['lower'][-1] - y['upper'][0]))
        # Contribution of panels not directly in contact with flow above the hinge
        Cm += (1./2*c**2)*(Cp['upper'][-1])*(((x['upper'][-1] + x_ref)/2. - \
              x_ref) *(x['upper'][-1] - x_ref) + ((y['upper'][-1] + y_ref)/2. - \
              y_ref)*(y['upper'][-1] - y_ref))
        # Contribution of panels not directly in contact with flow below the hinge
        Cm += (1./2*c**2)*(Cp['lower'][0])*(((x['lower'][0] + x_ref)/2. - x_ref) *(x_ref - x['lower'][0]) + \
              ((y['lower'][0] + y_ref)/2. - y_ref) *(y_ref - y['lower'][0]))
    return Cm
#==============================================================================
# Functions Intended for use with FInite ELement Methods
#==============================================================================
def force_shell(Data, chord, half_span, height, Velocity, thickness=0,
                txt=False):
    # Height is in feet
    # If the Shell is an extrude, it needs to take in consideration
    # that there is a skin thickness outwards of the outer mold.
    # If the Shell is na planar, there is no need for such a
    # consideration
    Air_properties = air_properties(height, unit='feet')
    atm_pressure = Air_properties['Atmospheric Pressure']
    air_density = Air_properties['Density']
    if thickness == 0:
        Data['Force'] = map(lambda Cp:(Cp*0.5*air_density * Velocity**2 +
                            atm_pressure) * chord*half_span, Data['Cp'])
        Data['x'] = map(lambda x: (chord)*x, Data['x'])
        Data['y'] = map(lambda x: (chord)*x, Data['y'])
    else:
        Data['Force'] = map(lambda Cp:(Cp*0.5*air_density * Velocity**2 +
                            atm_pressure) * chord*half_span, Data['Cp'])
        Data['x'] = map(lambda x: (chord - 2.*thickness) * x + thickness,
                        Data['x'])
        Data['y'] = map(lambda x: (chord - 2.*thickness) * x, Data['y'])
    Data['z'] = [0] * len(Data['x'])

    PressureDistribution = zip(Data['x'], Data['y'], Data['z'], Data['Force'])
#    elliptical_distribution=np.sqrt(1.-(Data['z']/half_span)**2)
#    if txt==True:
#        DataFile = open('Force_shell.txt','w')
#        DataFile.close()
#        for j in range(N):
#            for i in range(len(Data['x'])):
#                    DataFile = open('Force_shell.txt','a')
#                    DataFile.write('%f\t%f\t%f\t%f\n' % (
#                        Data['x'][i],
#                        Data['y'][i],
#                        Data['z'][j],
#                        elliptical_distribution[j]*Data['Force'][i]))
#                    DataFile.close()
#        return 0
#    else:
#        PressureDistribution=()
#        for j in range(N):
#            for i in range(len(Data['x'])):
#                    PressureDistribution=PressureDistribution+((Data['x'][i],
#                        Data['y'][i],Data['z'][j],
#                        elliptical_distribution[j]*Data['Pressure'][i]),)
    return PressureDistribution

def pressure_shell(Data, half_span, chord = 'MAX', air_density = 0, Velocity = 0,
                   N = 10, thickness = 0, txt=False, llt_distribution=False,
                   distribution='Uniform', amplifier = 1):
    """Converts pressure coefficient data, usually 2D, into a 3D presurre field
       that Abaqus understands. Can be used for shells (considers thicknesses),
       but also for any surface. Can do Lifting Line Theory (LLT), Elliptical,
       and Uniform distributions.

       If chord='MAX', the maximum value for vector 'x' is used as chord. If
       data in non-dimensional, use a numerical value.

       If txt==True, an output textfile is generated."""

    # If data is in the form of pressure coefficients, convert to pressure
    if 'Cp' in Data.keys():
        Data['Pressure'] = map(lambda Cp: Cp*0.5*air_density* Velocity**2 *chord,
                                Data['Cp'])
    if chord == 'MAX':
        chord = max(Data['x'])

    Data['x'] = map(lambda x: (chord - 2.*thickness)*x + thickness, Data['x'])
    Data['y'] = map(lambda x: (chord - 2.*thickness)*x, Data['y'])
    DataFile = open('Pressure_shell.txt', 'w')
    DataFile.close()
    if distribution == 'Elliptical':
        Data['z'] = np.linspace(0, half_span, N)
        distribution = amplifier*np.sqrt(1. - (Data['z']/half_span)**2)
    elif distribution == 'LLT':
        Data['z'] = np.linspace(0, half_span, N)
        distribution = amplifier*llt_distribution
    elif distribution == 'Uniform':
        Data['z'] = np.linspace(0, half_span, N)
        distribution = amplifier*np.ones((N,1))
    if txt == True:
        for j in range(N):
            for i in range(len(Data['x'])):
                DataFile = open('Pressure_shell.txt','a')
                DataFile.write('%f\t%f\t%f\t%f\n' % (
                    Data['x'][i],
                    Data['y'][i],
                    Data['z'][j],
                    distribution[j]*Data['Pressure'][i]))
                DataFile.close()
        return 0
    else:
        PressureDistribution = ()
        for j in range(N):
            for i in range(len(Data['x'])):
                    PressureDistribution = PressureDistribution + (
                                             (Data['x'][i], Data['y'][i],
                                              Data['z'][j], distribution[j]*
                                              Data['Pressure'][i]), )
        return PressureDistribution

def pressure_shell_2D(Data, chord, thickness, half_span, height, Velocity, N,
                      txt=False):
    """Calculate pressure field for a 2D Shell."""
    Air_properties = air_properties(height, unit='feet')
    air_density = Air_properties['Density']

    Data['Pressure'] = map(lambda Cp: Cp*0.5*air_density* Velocity**2 *chord,
                            Data['Cp'])
    Data['x'] = map(lambda x: (chord - 2.*thickness)*x + thickness, Data['x'])
    Data['y'] = map(lambda x: (chord - 2.*thickness)*x, Data['y'])
    DataFile = open('Pressure_shell.txt', 'w')
    DataFile.close()
    Data['z'] = np.linspace(0, half_span, N)
    if txt == True:
        for j in range(N):
            for i in range(len(Data['x'])):
                    DataFile = open('Pressure_shell.txt', 'a')
                    DataFile.write('%f\t%f\t%f\t%f\n' % (
                        Data['x'][i],
                        Data['y'][i],
                        Data['z'][j],
                        Data['Pressure'][i]))
                    DataFile.close()
        return 0
    else:
        PressureDistribution = ()
        for j in range(N):
            for i in range(len(Data['x'])):
                    PressureDistribution = PressureDistribution + (
                                            (Data['x'][i], Data['y'][i],
                                             Data['z'][j],
                                             Data['Pressure'][i]), )
        return PressureDistribution

def air_properties(height, unit='feet'):
    """ Function to calculate air properties for a given height (m or ft).

    Sources:
      - http://en.wikipedia.org/wiki/Density_of_air#Altitude
      - http://aerojet.engr.ucdavis.edu/fluenthelp/html/ug/node337.htm

    Created on Thu May 15 14:59:43 2014
    @author: Pedro Leal
    """
    # height is in m
    if unit == 'feet':
        height = 0.3048*height
    elif unit != 'meter':
        raise Exception('air_properties can onlu understand feet and meters')

    #==================================================================
    # Constants
    #==================================================================
    # Sea level standard atmospheric pressure
    P0 = 101325. # Pa
    # Sealevel standard atmospheric temperature
    T0 = 288.15 # K
    # Earth-surface gravitational acceleration
    g = 8.80655 # m/s2
    # Temperature lapse rate, 0.0065 K/m
    L = 0.0065 # K/m
    # Ideal (Universal) gas constant
    R = 8.31447 # J/(mol K)
    # Molar mass of dry air
    M = 0.0289644 #kg/mol
    # Specific R for air
    R_air = R/M
    # Sutherland's law coefficients
    C1 = 1.458e-6 #kg/m.s.sqrt(K)
    C2 = 110.4 #K

    #==================================================================
    # Temperature
    #==================================================================
    #Temperature at altitude h meters above sea level is approximated
    # by the following formula (only valid inside the troposphere):
    T = T0 - L*height

    #==================================================================
    # Pressure
    #==================================================================
    P = P0 * (1. - L*height/T0)**(g*M/(R*L))

    #==================================================================
    # Density
    #==================================================================
    density = P*M / (R*T)

    #==================================================================
    # Dynamic Viscosity (Sutherland equation with two constants)
    #==================================================================
    dyn_viscosity = (C1 * T**(3./2)) / (T+C2)

    return {'Density': density, 'Dynamic Viscosity': dyn_viscosity,
            'Atmospheric Temperature': T, 'R air': R_air,
            'Atmospheric Pressure': P}

def Reynolds(height, V, c):
    """Simple function to calculate Reynolds for a given height.

    @author: Pedro Leal
    Created in Jul 17 2015
    """

    Air_Data = air_properties(height, unit='feet')
    rho = Air_Data['Density']
    L = c
    nu = Air_Data['Dynamic Viscosity']
    return rho*V*L/nu

