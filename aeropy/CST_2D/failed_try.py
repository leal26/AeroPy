import math
import numpy as np
import matplotlib.pyplot as plt

def CST(x, c, deltasz=None, Au=None, Al=None, N1=0.5, N2=1., thetas = [0,0]):
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
        S=K(r,n)*(psi**r)*(1.-psi)**(n-r)
        return S
    
    # Class Function    
    def C(N1,N2,psi):
        C=((psi)**N1)*((1.-psi)**N2)
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
    C=C(N1,N2,psi);
    
    #==========================================================================
    #                   Defining the working surfaces
    #==========================================================================
    deltaz={}
    eta={}
    y={}
    x_r = {}
    y_r = {}
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
            n=len(A[surface])-1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Shape Function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
            Shape[surface]=0    
            for i in range(len(A[surface])):
#                print A
#                print S(i,n,psi)
                if surface=='l':
                    # if i==0:
                        # Shape[surface]-=A[surface][i]*S(i,n,psi)*math.cos(theta)
                    # else:
                        Shape[surface]-=A[surface][i]*S(i,n,psi)
                else:
                    # if i==0:
                        # Shape[surface]+=A[surface][i]*S(i,n,psi)*math.cos(theta)
                    # else:
                        Shape[surface]+=A[surface][i]*S(i,n,psi)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Airfoil Shape (eta=z/c)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
            # Airfoil Shape (eta=z/c)
            if surface=='l':
                eta[surface]= C*Shape[surface]-psi*deltaz[surface]/c #- \
                             #(1-psi)**n*(y_LE*math.cos(theta) + \
                             #x_LE*math.sin(theta));
            else:
                eta[surface]= C*Shape[surface]+psi*deltaz[surface]/c #+ \
                             #(1-psi)**n*(y_LE*math.cos(theta) + \
                             #x_LE*math.sin(theta));  
            # Giving back the dimensions
            y[surface]=c*eta[surface]
            
            # multiplying by angle function
            # theta_distribution = np.zeros(len(psi))
            # for i in range(len(thetas)):
            theta_distribution = thetas[0]*(1-psi)**n
            x_r[surface] = np.cos(theta_distribution)*x - np.sin(theta_distribution)*y[surface]
            y_r[surface] = np.sin(theta_distribution)*x + np.cos(theta_distribution)*y[surface]
            
    # redefining x and y for output
    x = x_r
    y = y_r
    if Al and Au:     
        return x,y
    elif Au:
        return x['u'],y['u']
    else:
        return x['l'],y['l']

def calculate_dependent_shape_coefficients(Au_C_1_to_n,
                                           psi_spars, Au_P, Al_P, deltaz, c_P,
                                           morphing = 'backwards', theta_C=0, theta_P=0):
    """Calculate  dependent shape coefficients for children configuration for a 4 order
    Bernstein polynomial and return the children upper, lower shape 
    coefficients, children chord and spar thicknesses. _P denotes parent parameters"""
    def calculate_AC_u0(AC_u0):
        Au_C = [AC_u0] + Au_C_1_to_n
        c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz)
        return np.sqrt(c_P/c_C)*Au_P[0]
    
    # Bersntein Polynomial
    def K(r,n):
        K=math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K
    # Bernstein Polynomial order
    n = len(Au_C_1_to_n)

    # Find upper shape coefficient though iterative method since Au_0 is unknown
    # via fixed point iteration
    #AC_u0 = optimize.fixed_point(calculate_AC_u0, Au_P[0])
    #print AC_u0
    error = 9999
    AC_u0 = Au_P[0]
    while error > 1e-9:
        before = AC_u0
        AC_u0 = calculate_AC_u0(AC_u0)
        error = abs(AC_u0 - before)

    # Because the output is an array, need the extra [0]      
    Au_C = [AC_u0] + Au_C_1_to_n
    
    # Now that AC_u0 is known we can calculate the actual chord and AC_l0
    c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz/c_P)
    AC_l0 = np.sqrt(c_P/c_C)*Al_P[0]
    # print '0 lower shape coefficient: ',AC_l0
    # Calculate thicknessed and tensor B for the constraint linear system problem
    spar_thicknesses = []
    A0 = AC_u0 + AC_l0
    
    if morphing == 'backwards':
        b_list = np.zeros((n,1))
        for j in range(len(psi_spars)):
            psi_j = psi_spars[j]
            #Calculate the spar thickness in meters from parent, afterwards, need to
            #adimensionalize for the goal airfoil by dividing by c_goal
            t_j = calculate_spar_distance(psi_spars[j], Au_C, Au_P, Al_P, deltaz, c_P)
            
            spar_thicknesses.append(t_j)
            b_list[j] = (t_j/c_C - psi_j*deltaz/c_C)/((psi_j**0.5)*(1-psi_j)) - A0*(1-psi_j)**n

        B = np.zeros((n,n))
        #j is the row dimension and i the column dimension in this case
        for j in range(n):
            for i in range(n):
                #Because in Python counting starts at 0, need to add 1 to be
                #coherent for equations
                r = i +1
                B[j][i] = K(r,n)*(psi_spars[j]**r)*(1-psi_spars[j])**(n-r)
        
        A_bar = np.dot(inv(B), b_list)

        Al_C = [AC_l0]
        for i in range(len(A_bar)):
            Al_C.append(A_bar[i][0] - Au_C[i+1]) #extra [0] is necessary because of array

    elif morphing == 'forwards':
        f = np.zeros((n,1))
        # psi/xi coordinates for lower surface of the children configuration
        psi_lower_children = []
        xi_lower_children = []
        xi_upper_children = []

        c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz)
        # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
        psi_upper_children = []
        for j in range(len(psi_spars)):
            psi_upper_children.append(calculate_psi_goal(psi_spars[j], Au_P, Au_C, deltaz,
                                   c_P, c_C))
        # Calculate xi for upper children. Do not care about lower so just gave it random shape coefficients
        xi_upper_children = CST(psi_upper_children, 1., deltasz= [deltaz/2./c_C, deltaz/2./c_C],  Al= Au_C, Au =Au_C, thetas=[theta_C])
        xi_upper_children = xi_upper_children['u']

        # print xi_upper_children
        
        #Debugging section
        x = np.linspace(0,1)
        y = CST(x, 1., deltasz= [deltaz/2./c_C, deltaz/2./c_C],  Al= Au_C, Au =Au_C, theta=[theta_C])
        # plt.plot(x,y['u'])
        # plt.scatter(psi_upper_children, xi_upper_children)
        # plt.grid()
        # plt.show()
        # BREAK
        for j in range(len(psi_spars)):
            xi_parent = CST(psi_spars, 1., deltasz= [deltaz/2./c_P, deltaz/2./c_P],  Al= Al_P, Au =Au_P, thetas=[theta_P])
            delta_j_P = xi_parent['u'][j]-xi_parent['l'][j]
            t_j = c_P*(delta_j_P)
            # Claculate orientation for children
            s_j = calculate_spar_direction(psi_spars[j], Au_P, Au_C, deltaz, c_C)
            psi_l_j = psi_upper_children[j]-delta_j_P/c_C*s_j[0]
            xi_l_j = xi_upper_children[j]-delta_j_P/c_C*s_j[1]

            spar_thicknesses.append(t_j)
            psi_lower_children.append(psi_l_j)
            xi_lower_children.append(xi_l_j)

            f[j] = (2*xi_l_j + psi_l_j*deltaz/c_C)/(2*(psi_l_j**0.5)*(psi_l_j-1))  - AC_l0*(1-psi_l_j)**n

        F = np.zeros((n,n))
        #j is the row dimension and i the column dimension in this case
        for j in range(n):
            for i in range(n):
                #Because in Python counting starts at 0, need to add 1 to be
                #coherent for equations
                r = i +1
                F[j][i] = K(r,n)*(psi_lower_children[j]**r)*(1-psi_lower_children[j])**(n-r)
        # print F
        # print f
        A_lower = np.dot(inv(F), f)

        Al_C = [AC_l0]
        for i in range(len(A_lower)):
            Al_C.append(A_lower[i][0]) #extra [0] is necessary because of array
    return Au_C, Al_C, c_C, spar_thicknesses

if __name__ == '__main__':        
    
    chord = 1
    Au = [0.2,0.1,0.1]
    Al = [0.2,0.1,0.1]
    deltaz = 0
    thetas = [0,0,.0]
    
    plt.figure()
    colors = ['b','g','r','k','c',]
    theta0_list = np.linspace(0,3.1415,3)
    #theta0_list = [0.75]
    for i in range(len(theta0_list)):
        thetas[0] = theta0_list[i]
        x = np.linspace(0, chord, 500)
        x,y = CST(x, chord, [deltaz/2., deltaz/2.], Au, Al, thetas = thetas)


        plt.plot(x['u'], y['u'], colors[i])
        plt.plot(x['l'], y['l'], colors[i], label=r'$\theta=$%.2f' % theta0_list[i] )
    plt.axis('equal')
    plt.scatter([0],[0],c='k')
    plt.grid()
    # plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    BREAK
    #==============================================================================
    # Structurally Consistent CST
    #==============================================================================
    # Parameter
    c_P = 1.                  #m
    deltaz = 0.*c_P    #m
    
    # Avian wing, order 5
    # Au_P = [0.23993240191629417, 0.34468227138908186, 0.18125405377549103, 
            # 0.35371349126072665, 0.2440815012119143, 0.25724974995738387]
    # Al_P = [0.18889012559339036, -0.24686758992053115, 0.077569769493868401,
            # -0.547827192265256, -0.0047342206759065641, -0.23994805474814629]
    # NACA0012
    Au_P =  [0.10887, 0.1187, 0.07843, 0.12084, 0.07919, 0.09840]
    Al_P =  [0.11117, 0.1000, 0.1239, 0.06334, 0.11539, 0.10400]  
    n = len(Au_P) - 1
    theta_P = 0
    if inverted:
        temp = Au_P
        Au_P = list(-np.array(Al_P))
        Al_P = list(-np.array(temp))
    # Shape coefficients for upper surface of cruise airfoil
    # Small
    AC_u1 = 0.1487            #Adimensional
    AC_u2 = 0.10843          #Adimensional
    AC_u3 = 0.15084                #Adimensional
    AC_u4 = 0.10919              #Adimensional
    AC_u5 = 0.12840             #Adimensional
    theta_C = 0.5
    #Spar position for cruise (adiminesional because the chord will still be calculated)
    psi_spar1 = 0.2           #Adimensional
    psi_spar2 = 0.3           #Adimensional
    psi_spar3 = 0.5                                         #Adimensional
    psi_spar4 = 0.7                                         #Adimensional
    psi_spar5 = 0.9                                         #Adimensional
    psi_spars = [psi_spar1, psi_spar2, psi_spar3, psi_spar4, psi_spar5]

    #==============================================================================
    # Calculate dependent coefficients
    #==============================================================================
    Au_C, Al_C, c_C, spar_thicknesses = calculate_dependent_shape_coefficients(
                                                        AC_u1, AC_u2, AC_u3, AC_u4, AC_u5,
                                                        psi_spars, Au_P, Al_P,
                                                        deltaz, c_P, morphing=morphing_direction,
                                                        theta_C=theta_C, theta_P=theta_P)
    
    #==============================================================================
    #  Plot results
    #==============================================================================
    np.set_printoptions(precision=20)
    # Print shape for children
    x = np.linspace(0, c_C, 100000)
    y = CST(x, c_C, deltasz= [deltaz/2., deltaz/2.],  Al= Al_C, Au =Au_C, thetas=[theta_C])

    plt.plot(x, y['u'], 'b', label = 'Children', lw=2)
    plt.plot(x, y['l'], 'b', label = None, lw=2)


    # Print shape for parent
    x = np.linspace(0, c_P, 100000)
    y = CST(x, c_P, deltasz= [deltaz/2., deltaz/2.],  Al= Al_P, Au =Au_P, thetas=[theta_P])
    plt.plot(x, y['u'], 'r--', label='Parent', lw=2)
    plt.plot(x, y['l'], 'r--', label = None, lw=2)
    
    if morphing_direction == 'forwards':
        psi_flats = []
        intersections_x_children = [0]
        intersections_y_children = [0]
        intersections_x_parent = [0]
        intersections_y_parent = [0]
        for j in range(len(psi_spars)):
            psi_parent_j = psi_spars[j]
            # Calculate psi at landing
            # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
            psi_children_j = calculate_psi_goal(psi_parent_j, Au_P, Au_C, deltaz, c_P, c_C)
            x_children_j = psi_children_j*c_C
            
            # Calculate xi at landing
            temp = CST(x_children_j, c_C, [deltaz/2., deltaz/2.], Al= Al_C, Au =Au_C, thetas=[theta_C])
            y_children_j = temp['u']
        
            s = calculate_spar_direction(psi_spars[j], Au_P, Au_C, deltaz, c_C)
            
            # Print spars for children
            if not inverted:
                plt.plot([x_children_j, x_children_j - spar_thicknesses[j]*s[0]],[y_children_j, y_children_j - spar_thicknesses[j]*s[1]], c = 'b', lw=2, label=None)
            else:
                plt.plot([x_children_j, x_children_j - spar_thicknesses[j]*s[0]],[-y_children_j, -y_children_j + spar_thicknesses[j]*s[1]], c = 'b', lw=2, label=None)
            psi_flats.append(x_children_j - spar_thicknesses[j]*s[0])
            y = CST(np.array([psi_parent_j*c_P]), c_P, deltasz=[deltaz/2., deltaz/2.], Al= Al_P, Au =Au_P, thetas=[theta_P])
            
            intersections_x_children.append(x_children_j - spar_thicknesses[j]*s[0])
            intersections_y_children.append(y_children_j - spar_thicknesses[j]*s[1])
            
            # Print spars for parents
            if not inverted:
                plt.plot([psi_parent_j*c_P, psi_parent_j*c_P], [y['u'], y['u']-spar_thicknesses[j]], 'r--', lw=2, label = None)
            else:
                plt.plot([psi_parent_j*c_P, psi_parent_j*c_P], [-y['u'], -y['u']+spar_thicknesses[j]], 'r--', lw=2, label = None)

            intersections_x_parent.append(psi_parent_j*c_P)
            intersections_y_parent.append(y['u']-spar_thicknesses[j])