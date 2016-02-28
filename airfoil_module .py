# -*- coding: utf-8 -*-
"""
Created on Fri Oct 09 17:31:59 2015

Current functionalities:
- CST
- flaps
- rotate
@author: Pedro Leal
"""

import math
import numpy as np
    
def CST(x,c,deltasz=None,Au=None,Al=None):
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

        - x:list of points along the chord, from TE and the LE, or vice-versa.
          The code works both ways.
          
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
    def K(r,n):
        K=math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K
    
    # Shape Function   
    def S(r,n,psi):
        S=K(r,n)*(psi**r)*(1-psi)**(n-r)
        return S
    
    # Class Function    
    def C(N1,N2,psi):
        C=((psi)**N1)*((1-psi)**N2)
        return C
    
    if type(x)==list:
        x=np.array(x)
    # Adimensionalizing
    psi=x/c;    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Class Function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # The Coefficients for an airfoil with a rounded leading edge and a sharp
    # trailing edge are N1=0.5 and N2=1.0.
    N1=0.5;
    N2=1.0;
    C=C(N1,N2,psi);
    
    #==========================================================================
    #                   Defining the working surfaces
    #==========================================================================
    deltaz={}
    eta={}
    y={}
    Shape={}
    
    if Al and Au:
        deltaz['u']=deltasz[0]
        deltaz['l']=deltasz[1]
        
        if len(Au)!=len(Al):
            raise Exception("Au and Al need to have the same dimensions")
        elif len(deltasz)!=2:
            raise Exception("If both surfaces are being analyzed, two values for deltasz are needed")
        
    elif Au and not Al:
        if type(deltasz)==list:
            if len(deltaz['u'])!=1:
                raise Exception("If only one surface is being analyzed, one value for deltasz is needed")
            else:
                deltaz['u']=float(deltasz)
        else:
            deltaz['u']=deltasz
       
    elif Al and not Au:
        if type(deltasz)==list:
            if (deltaz['l'])!=1:
                raise Exception("If only one surface is being analyzed, one value for deltasz is needed")
            else:
                deltaz['l']=float(deltasz)
        else:
            deltaz['l']=deltasz
    else:
        raise Exception("Au or Al need to have at least one value")
    A={'u':Au,'l':Al}
    for surface in ['u','l']:
        if A[surface]:
            
            if type(A[surface])==int or type(A[surface])==float:
                A[surface]=[A[surface]]
            # the degree of the Bernstein polynomial is given by the number of
            # coefficients
            n=len(A[surface])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Shape Function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
            Shape[surface]=0    
            for i in range(len(A[surface])):
#                print A
#                print S(i,n,psi)
                if surface=='l':
                    Shape[surface]-=A[surface][i]*S(i,n,psi)
                else:
                    Shape[surface]+=A[surface][i]*S(i,n,psi)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Airfoil Shape (eta=z/c)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
            # Airfoil Shape (eta=z/c)
            if surface=='l':
                eta[surface]= C*Shape[surface]-psi*deltaz[surface]/c;
            else:
                eta[surface]= C*Shape[surface]+psi*deltaz[surface]/c;  
            # Giving back the dimensions
            y[surface]=c*eta[surface]
    if Al and Au:     
        return y
    elif Au:
        return y['u']
    else:
        return y['l']

#==============================================================================
# The following function are related to the use of plain flaps
#==============================================================================
    
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
    
            #If equal just save and break it
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
        
    #Find y hinge
    upper_fwd = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    upper_aft = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    lower_fwd = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    lower_aft = {'x': 9e9, 'y': 9e9, 'error': 9e9}
    hinge = {'x': x_hinge, 'y': 9e9}
    
    upper_fwd, upper_aft = find_closest(upper, x_hinge, upper_fwd, upper_aft)
    lower_fwd, lower_aft = find_closest(lower, x_hinge, lower_fwd, lower_aft) 
    
    #Interpolate
    hinge['y_upper'] =  upper_fwd['y'] + (
                        hinge['x'] - upper_fwd['x'])*(
                        upper_aft['y'] - upper_fwd['y'])/(
                        upper_aft['x'] - upper_fwd['x'])
    
    hinge['y_lower'] =  lower_fwd['y'] + (
                        hinge['x'] - lower_fwd['x'])*(
                        lower_aft['y'] - lower_fwd['y'])/(
                        lower_aft['x'] - lower_fwd['x'])
    hinge['y']  = (hinge['y_upper'] + hinge['y_lower'])/2.

    return hinge


def find_flap(data, hinge):
    """Create the static airfoil and flap dictionaries containing the outer
    mold coordinates of both.
    
    :param data: dictionary with x and y cooridnates of the whole outer mold
    
    :param hinge: dictionary with x and y coordinates of hinge. Can be found
                  via find_hinge."""
    #Generate empty dictionaries which will store final data.
    flap_data = {'x':[], 'y':[]}
    static_data = {'x':[], 'y':[]}
    
    type = None
    
    for i in range(len(data['x'])):
        xi = data['x'][i]
        yi = data['y'][i]
    #Because of the way that xfoil works, the upper list will always
    #begin from the trailing edge    
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
                
            flap_data['x'].append(xi)
            flap_data['y'].append(yi)
            
    #Because of the way that xfoil works, the lower list will always
    #begin from the  leading edge  
        else:
            if type == None:
                type = 'lower'
            elif type == 'upper':
                flap_data['x'].append(hinge['x'])
                flap_data['y'].append(hinge['y_' + type])
                static_data['x'].append(hinge['x'])
                static_data['y'].append(hinge['y_' + type])
                # this will make the at the hinge to be included only once
                type = 'lower'
                
            static_data['x'].append(xi)
            static_data['y'].append(yi)
    return static_data, flap_data

def rotate(upper, lower, origin, theta):
    """
    :param upper: dictionary with keys x and y, each a list
    
    :param lower: dictionary with heys x and y, each a list
    
    :param origin: dictionary with keys x and y, each a float
    
    :param theta: float representing angle in degrees clock-wise
    """
    output = []
    
    #For trigonometric relations in numpy, theta must be in radians
    theta = theta * np.pi/180
    # Rotation transformation Matrix
    T = [[np.cos(theta), np.sin(theta)],
        [-np.sin(theta), np.cos(theta)]]
        
    for coordinates in [upper, lower]:
        rotated_coordinates = {'x':[], 'y':[]}
        for i in range(len(coordinates['x'])):
            # The rotation must take place with the center of rotation
            # as the origin
            cx = coordinates['x'][i] - origin['x']
            cy = coordinates['y'][i] - origin['y']
            # Rotate
            rot_x = T[0][0]*cx + T[0][1]*cy
            rot_y = T[1][0]*cx + T[1][1]*cy
            # Store and add back the values of the origin
            rotated_coordinates['x'].append(rot_x + origin['x'])
            rotated_coordinates['y'].append(rot_y + origin['y'])
        output.append(rotated_coordinates)
    return output[0], output[1]

def clean(upper_static, upper_flap, lower_static, lower_flap, hinge, 
          deflection, N):
    """
    Function to remove intersecting lines and to round transition
    on upper surface.
    
    :param N: number of points on smooth transition
    """
    #The lower surface has a non-smooth trnaisiton so we only have to
    # clean the intersecting lines
    intersection= intersect_curves(lower_static['x'], lower_static['y'],
                                   lower_flap['x'], lower_flap['y'],
                                   input_type = 'list')
    print intersection
    intersection_x = intersection[0][0]
    intersection_y = intersection[1][0]
          
    for i in range(len(lower_flap['x'])):
        xi = lower_flap['x'][i]
        # From the plots, every point of the flap before the intersection
        # needs to be elimnated.
        if xi < intersection_x:
            lower_flap['x'][i] = None
            lower_flap['y'][i] = None
            
    for i in range(len(lower_static['x'])):
        xi = lower_static['x'][i]
        # From the plots, every point of the flap before the intersection
        # needs to be elimnated.
        if xi > intersection_x:
            lower_static['x'][i] = None
            lower_static['y'][i] = None
            
    #Eliminatting the None vectors
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
    # Need to create points connecting upper part in circle (part of plain flap
    # exposed during flight)
    for i in range(N):
        theta = ((N-1-i)/(N-1.)) * deflection * np.pi/180.
        upper_smooth['x'].append(hinge['x'] + R*np.sin(theta))
        upper_smooth['y'].append(hinge['y'] + R*np.cos(theta))
    # Assembling all together
    
    modified_airfoil = {'x':[], 'y':[]}
    for key in ['x','y']:
        modified_airfoil[key] = upper_flap[key] + upper_smooth[key] + \
                                upper_static[key] + lower_static[key] + \
                                lower_flap[key]
    return modified_airfoil

def intersect_curves(x1, y1, x2, y2, input_type = 'list'):
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
    #convert inputs to numpy arrays
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
    mix1 = np.tile(np.amin(xs1, axis=0), (N2-1,1))
    max1 = np.tile(np.amax(xs1, axis=0), (N2-1,1))
    miy1 = np.tile(np.amin(ys1, axis=0), (N2-1,1))
    may1 = np.tile(np.amax(ys1, axis=0), (N2-1,1))
    mix2 = np.transpose(np.tile(np.amin(xs2, axis=0), (N1-1,1)))
    max2 = np.transpose(np.tile(np.amax(xs2, axis=0), (N1-1,1)))
    miy2 = np.transpose(np.tile(np.amin(ys2, axis=0), (N1-1,1)))
    may2 = np.transpose(np.tile(np.amax(ys2, axis=0), (N1-1,1)))
    idx = np.where((mix2 <= max1) & (max2 >= mix1) & (miy2 <= may1) & (may2 >= miy1)) # overlapping segment combinations

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
        a = np.array([[xa[1] - xa[0], xb[0] - xb[1]], [ya[1] - ya[0], yb[0]- yb[1]]])
        b = np.array([xb[0] - xa[0], yb[0] - ya[0]])
        r, residuals, rank, s = np.linalg.lstsq(a, b)
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

    return (x0, y0)

def find_deflection(x_hinge, upper_cruise, lower_cruise, 
                    type = 'Match trailing edge', alpha=0, Reynolds = 0, 
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
    
    #==============================================================================
    #  Find Deflection      
    #==============================================================================
    if type == 'Match trailing edge':
        #Calculate deflection angle in radians
        deflection = np.arctan2(hinge['y'] - y_TE_objective, x_TE_objective - hinge['x']) - np.arctan2(hinge['y'] - y_TE_baseline, x_TE_baseline - hinge['x'])
        # Convet to degrees
        deflection = deflection*180./np.pi
    
        upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)
        
        flapped_airfoil = clean(upper_static, upper_rotated, lower_static, 
                                lower_rotated, hinge, deflection, N = 5)
    
        xf.create_input(airfoil, 'Type 3', airfoil = flapped_airfoil)
                                    
        xf.call(airfoil, alfas = alpha, Reynolds = Reynolds, iteration=100, GDES = True,
                output='Polar')
        filename = xf.file_name(airfoil, alpha, output='Polar')
        Data = xf.output_reader(filename, output='Polar')
        
        #Rename to make more convenient
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
        while current_Cl < Cl_objective and deflection < 90.:# Convet to degrees
            previous_Cl = current_Cl
            previous_Cd = current_Cd
            previous_deflection = deflection
            history.append([previous_deflection, previous_Cl, previous_Cd])
            deflection = deflection + step
            print deflection
            upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)
            
            flapped_airfoil = clean(upper_static, upper_rotated, lower_static, 
                                    lower_rotated, hinge, deflection, N = 5)
        
            xf.create_input(airfoil, 'Type 3', airfoil = flapped_airfoil)
                                        
            xf.call(airfoil, alfas = alpha, Reynolds = Reynolds, iteration=100, GDES = True,
                    output='Polar')
            filename = xf.file_name(airfoil, alpha, output='Polar')
            Data = xf.output_reader(filename, output='Polar')
            
            #Rename to make more convenient
            current_Cd = Data['CD'][0]
            current_Cl = Data['CL'][0]

        for i in range(1,6):
            # Fine search
            step = 0.5**i
            current_Cl = previous_Cl
            while current_Cl < Cl_objective and deflection < 90.:# Convet to degrees
                previous_Cl = current_Cl
                previous_Cd = current_Cd
                previous_deflection = deflection
                history.append([previous_deflection, previous_Cl, previous_Cd])
                deflection = deflection + step
            
                upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)
                previous_flapped_airfoil = flapped_airfoil
                flapped_airfoil = clean(upper_static, upper_rotated, lower_static, 
                                        lower_rotated, hinge, deflection, N = 5)
            
                xf.create_input(airfoil, 'Type 3', airfoil = flapped_airfoil)
                                            
                xf.call(airfoil, alfas = alpha, Reynolds = Reynolds, iteration=100, GDES = True,
                        output='Polar')
                filename = xf.file_name(airfoil, alpha, output='Polar')
                Data = xf.output_reader(filename, output='Polar')
                
                #Rename to make more convenient
                current_Cd = Data['CD'][0]
                current_Cl = Data['CL'][0]
            print current_Cl, deflection

    elif type == 'Same Cd':
        current_Cl = 0
        current_Cd = 0
        previous_Cl = 0
        deflection = 2.
        history = []
        # Rough search
        step = 2.
        while current_Cd < Cd_objective and deflection < 90.:# Convet to degrees
            previous_Cl = current_Cl
            previous_Cd = current_Cd
            previous_deflection = deflection
            history.append([previous_deflection, previous_Cl, previous_Cd])
            deflection = deflection + step
            print deflection
            upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)
            
            flapped_airfoil = clean(upper_static, upper_rotated, lower_static, 
                                    lower_rotated, hinge, deflection, N = 5)
        
            xf.create_input(airfoil, 'Type 3', airfoil = flapped_airfoil)
                                        
            xf.call(airfoil, alfas = alpha, Reynolds = Reynolds, iteration=100, GDES = True,
                    output='Polar')
            filename = xf.file_name(airfoil, alpha, output='Polar')
            Data = xf.output_reader(filename, output='Polar')
            
            #Rename to make more convenient
            current_Cd = Data['CD'][0]
            current_Cl = Data['CL'][0]

        for i in range(1,6):
            # Fine search
            step = 0.5**i
            current_Cd = previous_Cd
            while current_Cd < Cd_objective or deflection < 90.:# Convet to degrees
                previous_Cl = current_Cl
                previous_Cd = current_Cd
                previous_deflection = deflection
                history.append([previous_deflection, previous_Cl, previous_Cd])
                deflection = deflection + step
            
                upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)
                previous_flapped_airfoil = flapped_airfoil
                flapped_airfoil = clean(upper_static, upper_rotated, lower_static, 
                                        lower_rotated, hinge, deflection, N = 5)
            
                xf.create_input(airfoil, 'Type 3', airfoil = flapped_airfoil)
                                            
                xf.call(airfoil, alfas = alpha, Reynolds = Reynolds, iteration=1000, GDES = True,
                        output='Polar')
                filename = xf.file_name(airfoil, alpha, output='Polar')
                Data = xf.output_reader(filename, output='Polar')
                
                #Rename to make more convenient
                print deflection
                current_Cd = Data['CD'][0]
                current_Cl = Data['CL'][0]
            print current_Cd, deflection
            print history
            print Cd_objective
        deflection = previous_deflection
        Cd = previous_Cd
        Cl = previous_Cl
        flapped_airfoil = previous_flapped_airfoil
        
    elif type == 'Range of deflection':
        current_Cl = 0
        current_Cd = 0
        previous_Cl = 0
        deflection = init_deflection
        history = {'CD':[], 'CL':[], 'deflection':[]}
        # Rough search
        while deflection < max_deflection:# Convet to degrees
            previous_Cl = current_Cl
            previous_Cd = current_Cd
            previous_deflection = deflection

            deflection = deflection + step
            print deflection
            upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)
            
            flapped_airfoil = clean(upper_static, upper_rotated, lower_static, 
                                    lower_rotated, hinge, deflection, N = 5)
        
            xf.create_input(airfoil, 'Type 3', airfoil = flapped_airfoil)
                                        
            xf.call(airfoil, alfas = alpha, Reynolds = Reynolds, iteration=300, GDES = True,
                    output='Polar')
            filename = xf.file_name(airfoil, alpha, output='Polar')
            Data = xf.output_reader(filename, output='Polar')
            
            #Rename to make more convenient
            print deflection
            current_Cd = Data['CD'][0]
            current_Cl = Data['CL'][0]
            if deflection <= max_deflection:
                history['CD'].append(current_Cd)
                history['CL'].append(current_Cl)
                history['deflection'].append(deflection)
        return history
    return deflection, Cd, Cl, flapped_airfoil

if __name__ == '__main__':
    import matplotlib.pyplot as plt    
    
    import xfoil_module as xf

    airfoil = "naca0012"
    xf.call(airfoil, output='Coordinates')
    filename = xf.file_name(airfoil, output='Coordinates')
    Data = xf.output_reader(filename, output='Coordinates', header=['x','y'])

    def separate_upper_lower(x,y):
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
        return x, y,
    
    x, y = separate_upper_lower(Data['x'],Data['y'])
    
    upper = {'x': x['upper'], 'y': y['upper']}
    lower = {'x': x['lower'], 'y': y['lower']}
    
 
    x_hinge = 0.75
    hinge = find_hinge(x_hinge, upper, lower)
    
    upper_static, upper_flap = find_flap(upper, hinge)
    lower_static, lower_flap = find_flap(lower, hinge)
#    plt.scatter(upper_rotated['x'], upper_rotated['y'])
#    plt.scatter(lower_static['x'], lower_static['y'])    
    deflection = 5.4 # Deflection at which the flap and morphed have the same drag
    
    upper_rotated, lower_rotated = rotate(upper_flap, lower_flap, hinge, deflection)
    
    plt.scatter(upper_rotated['x'], upper_rotated['y'])
    plt.scatter(lower_static['x'], lower_static['y'], c = 'r')
    plt.scatter(lower_rotated['x'], lower_rotated['y'])
    plt.scatter(upper_static['x'], upper_static['y'], c = 'r')
    flapped_airfoil = clean(upper_static, upper_rotated, lower_static, 
                            lower_rotated, hinge, deflection, N = 5)
    