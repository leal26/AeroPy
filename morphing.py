# -*- coding: utf-8 -*-
"""
Objective: create an airfoil with a leading edge restriction, same upper length
restriction, othogonal upper spars and constant thicknesses in four places

Created on Mon Oct 17 10:36:34 2016

@author: Pedro
"""
import math
import numpy as np
from numpy.linalg import inv

from airfoil_module import CST
from CST_module import *

# Just as quick trick, to make upper morph I just mirror the image in regards to x
inverted = False
# Defines if basckwards or forwards morphing
morphing_direction = 'forwards'
	
#==============================================================================
# Calculate dependent shape function parameters
#==============================================================================
def calculate_dependent_shape_coefficients(AC_u1, AC_u2, AC_u3, AC_u4, AC_u5,
                                           psi_spars, Au_P, Al_P, deltaz, c_P,
                                           morphing = 'backwards'):
    """Calculate  dependent shape coefficients for children configuration for a 4 order
    Bernstein polynomial and return the children upper, lower shape 
    coefficients, children chord and spar thicknesses. _P denotes parent parameters"""
    def calculate_AC_u0(AC_u0):
        Au_C = [AC_u0, AC_u1, AC_u2, AC_u3, AC_u4, AC_u5]
        c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz)
        return np.sqrt(c_P/c_C)*Au_P[0]
    
    # Bersntein Polynomial
    def K(r,n):
        K=math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K
    # Bernstein Polynomial order
    n = 5

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
    Au_C = [AC_u0, AC_u1, AC_u2, AC_u3, AC_u4, AC_u5]
    
    # Now that AC_u0 is known we can calculate the actual chord and AC_l0
    c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz/c_P)
    AC_l0 = np.sqrt(c_P/c_C)*Al_P[0]
    print '0 lower shape coefficient: ',AC_l0
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
        xi_upper_children = CST(psi_upper_children, 1., deltasz= [deltaz/2./c_C, deltaz/2./c_C],  Al= Au_C, Au =Au_C)
        xi_upper_children = xi_upper_children['u']

        print xi_upper_children
        
        #Debugging section
        x = np.linspace(0,1)
        y = CST(x, 1., deltasz= [deltaz/2./c_C, deltaz/2./c_C],  Al= Au_C, Au =Au_C)
        # plt.plot(x,y['u'])
        # plt.scatter(psi_upper_children, xi_upper_children)
        # plt.grid()
        # plt.show()
        # BREAK
        for j in range(len(psi_spars)):
            xi_parent = CST(psi_spars, 1., deltasz= [deltaz/2./c_P, deltaz/2./c_P],  Al= Al_P, Au =Au_P)
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
        print F
        print f
        A_lower = np.dot(inv(F), f)

        Al_C = [AC_l0]
        for i in range(len(A_lower)):
            Al_C.append(A_lower[i][0]) #extra [0] is necessary because of array
    return Au_C, Al_C, c_C, spar_thicknesses

def calculate_shape_coefficients_tracing(A0, tip_displacement, other_points, N1, N2):   
    """
    inputs:
        - tip_displacement: {'x': value, 'y': value}
        - other_points: {'x': value, 'y': value}
        - A0: float value for first shape coefficient. Usually related to a constraint.
    """
    # Bersntein Polynomial
    def K(r,n):
        K=math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K
 
    n = len(other_points['x'])
    chord = tip_displacement['y']
    
    Psi = np.array(other_points['y'])/chord
    Xi = np.array(other_points['x'])/chord
    
    EndThickness = tip_displacement['x']/chord
	
    T = np.zeros((n,n))
    t = np.zeros((n,1))
    for j in range(1,n+1):
        jj = j - 1
        for i in range(1,n+1):
            print i,j
            ii = i -1
            T[jj][ii] = K(i,n)* Psi[jj]**i * (1-Psi[jj])**(n-i)
        t[jj] = (Xi[jj] - Psi[jj]*EndThickness)/(Psi[jj]**N1*(1-Psi[jj])**N2) - A0*(1-Psi[jj])**n
    print T
    print t
    print (Xi[ii] - Psi[ii]*EndThickness)/(Psi[ii]*(1-Psi[ii])) - A0
    # Calculate the inverse
    A = np.dot(inv(T), t)
    A = [A0] + list(A.transpose()[0])
    print A
    return A
    
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    #testing = 'structurally_consistent'
    testing = 'tracing'
    
    if testing == 'tracing':
        N1 = 1.
        N2 = 1.
        tip_displacement = {'x': .1, 'y':1.}
        other_points = {'x': [0.01, -0.03, .05, 0.12], 'y':[0.1, 0.3, .5, 0.8]}
        A0 = -tip_displacement['x']
        print A0
        A = calculate_shape_coefficients_tracing(A0, tip_displacement, other_points, N1, N2)
        
        #plotting
        y = np.linspace(0, tip_displacement['y'], 100000)
        x = CST(y, tip_displacement['y'], deltasz= tip_displacement['x'],  Au = A, N1=N1, N2=N2)
        plt.plot(x,y)
        plt.scatter(other_points['x'] + [tip_displacement['x']], 
                    other_points['y'] + [tip_displacement['y']])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
        
    elif testing == 'structurally_consistent':
        
        #==============================================================================
        # Inputs
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
        
        if inverted:
            temp = Au_P
            Au_P = list(-np.array(Al_P))
            Al_P = list(-np.array(temp))
        # Shape coefficients for upper surface of cruise airfoil
        # AC_u1 = 0.25           #Adimensional
        # AC_u2 = 0.25          #Adimensional
        # AC_u3 = 0.25                #Adimensional
        # AC_u4 = 0.25             #Adimensional
        # AC_u5 = 0.25   	
        # Medium
        # AC_u1 = 0.2187            #Adimensional
        # AC_u2 = 0.17843          #Adimensional
        # AC_u3 = 0.22084                #Adimensional
        # AC_u4 = 0.17919              #Adimensional
        # AC_u5 = 0.19840             #Adimensional
        # Small
        AC_u1 = 0.1487            #Adimensional
        AC_u2 = 0.10843          #Adimensional
        AC_u3 = 0.15084                #Adimensional
        AC_u4 = 0.10919              #Adimensional
        AC_u5 = 0.12840             #Adimensional
        
        # AC_u1 = 0.34468227138908186                #Adimensional
        # AC_u2 = 0.18125405377549103                 #Adimensional
        # AC_u3 = 0.35371349126072665                #Adimensional
        # AC_u4 = 0.2440815012119143                 #Adimensional
        # AC_u5 = 0.25724974995738387                 #Adimensional
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
                                                            deltaz, c_P, morphing=morphing_direction)
        
        #==============================================================================
        #  Plot results
        #==============================================================================
        np.set_printoptions(precision=20)
        # Print shape for children
        x = np.linspace(0, c_C, 100000)
        y = CST(x, c_C, deltasz= [deltaz/2., deltaz/2.],  Al= Al_C, Au =Au_C)

        plt.plot(x, y['u'], 'b', label = 'Children', lw=2)
        plt.plot(x, y['l'], 'b', label = None, lw=2)


        # Print shape for parent
        x = np.linspace(0, c_P, 100000)
        y = CST(x, c_P, deltasz= [deltaz/2., deltaz/2.],  Al= Al_P, Au =Au_P)
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
                temp = CST(x_children_j, c_C, [deltaz/2., deltaz/2.], Al= Al_C, Au =Au_C)
                y_children_j = temp['u']
            
                s = calculate_spar_direction(psi_spars[j], Au_P, Au_C, deltaz, c_C)
                
                # Print spars for children
                if not inverted:
                    plt.plot([x_children_j, x_children_j - spar_thicknesses[j]*s[0]],[y_children_j, y_children_j - spar_thicknesses[j]*s[1]], c = 'b', lw=2, label=None)
                else:
                    plt.plot([x_children_j, x_children_j - spar_thicknesses[j]*s[0]],[-y_children_j, -y_children_j + spar_thicknesses[j]*s[1]], c = 'b', lw=2, label=None)
                psi_flats.append(x_children_j - spar_thicknesses[j]*s[0])
                y = CST(np.array([psi_parent_j*c_P]), c_P, deltasz=[deltaz/2., deltaz/2.], Al= Al_P, Au =Au_P)
                
                intersections_x_children.append(x_children_j - spar_thicknesses[j]*s[0])
                intersections_y_children.append(y_children_j - spar_thicknesses[j]*s[1])
                
                # Print spars for parents
                if not inverted:
                    plt.plot([psi_parent_j*c_P, psi_parent_j*c_P], [y['u'], y['u']-spar_thicknesses[j]], 'r--', lw=2, label = None)
                else:
                    plt.plot([psi_parent_j*c_P, psi_parent_j*c_P], [-y['u'], -y['u']+spar_thicknesses[j]], 'r--', lw=2, label = None)

                intersections_x_parent.append(psi_parent_j*c_P)
                intersections_y_parent.append(y['u']-spar_thicknesses[j])
        elif morphing_direction == 'backwards':
            # For backwards, goal is the parent and deformed is children
            for i in range(len(psi_spars)):
                psi_i = psi_spars[i]
                # Calculate psi at landing
                psi_goal_i = calculate_psi_goal(psi_i, Au_C, Au_P, deltaz, c_C, c_P)
                x_goal_i = psi_goal_i*c_P
                # Calculate xi at landing
                temp = CST(x_goal_i, c_P, [deltaz/2., deltaz/2.], Al= Al_P, Au =Au_P)
                y_goal_i = temp['u']

                #calculate spar direction
                s = calculate_spar_direction(psi_i, Au_C, Au_P, deltaz, c_P)

                plt.plot([x_goal_i, x_goal_i - spar_thicknesses[i]*s[0]],[y_goal_i, y_goal_i - spar_thicknesses[i]*s[1]], 'r--')

                y = CST(np.array([psi_i*c_C]), c_C, deltasz=[deltaz/2., deltaz/2.], Al= Al_C, Au =Au_C)

                plt.plot([psi_i*c_C, psi_i*c_C], [y['u'], y['u']-spar_thicknesses[i]], 'b', lw=2, label = None)

        plt.xlabel('$\psi^p$', fontsize = 14)
        plt.ylabel(r'$\xi^p$', fontsize = 14)
        plt.ylim([-0.06,0.17])
        plt.grid()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.legend(loc=1)
        plt.show()
        
        if morphing_direction == 'forwards':
            print c_C, c_P
            # Calculate initial lengths
            psi_list = [0.] + psi_spars + [c_P]
            print psi_list
            initial_lengths = []
            for i in range(len(psi_list)-1):
                initial_lengths.append(calculate_arc_length(psi_list[i], psi_list[i+1], Al_P, deltaz, c_P))
            # Calculate final lengths
            final_lengths = []
            psi_list = [0.] + psi_flats + [c_C] # In P configuration
            print psi_list
            for i in range(len(psi_list)-1):
                print psi_list[i]*c_P/c_C, psi_list[i+1]*c_P/c_C
                final_lengths.append(calculate_arc_length(psi_list[i]*c_P/c_C, psi_list[i+1]*c_P/c_C, Al_C, deltaz, c_C))
            # Calculate strains
            strains = []
            for i in range(len(final_lengths)):
                strains.append((final_lengths[i]-initial_lengths[i])/initial_lengths[i])
            
            for i in range(len(strains)):
                print 'Initial length: ' + str(initial_lengths[i]) + ', final length: ' + str(final_lengths[i]) + ', strains: ' + str(strains[i])
                
            intersections_x_children.append(c_C)
            intersections_y_children.append(0)
            intersections_x_parent.append(c_P)
            intersections_y_parent.append(0)		
            # Wire lengths
            for i in range(len(intersections_x_children)-1):
                length_parent = math.sqrt((intersections_x_parent[i]-intersections_x_parent[i+1])**2+
                                          (intersections_y_parent[i]-intersections_y_parent[i+1])**2)
                length_children = math.sqrt((intersections_x_children[i]-intersections_x_children[i+1])**2+
                                            (intersections_y_children[i]-intersections_y_children[i+1])**2)
                print (length_children-length_parent)/length_parent