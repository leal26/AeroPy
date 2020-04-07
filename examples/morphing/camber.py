import numpy as np
import matplotlib.pyplot as plt
from aeropy.geometry.airfoil import CST
from aeropy.morphing.camber_2D import *
# testing = 'structurally_consistent'
testing = 'structurally_consistent'
inverted = False
morphing_direction = 'forwards'
if testing == 'tracing':
    N1 = 1.
    N2 = 1.
    tip_displacement = {'x': 1., 'y':.5}
    other_points = {'x': [0.7], 'y':[0.25]}
    A0 = -tip_displacement['x']/tip_displacement['y']
    
    # Check if y values are smaller than tip y
    for y_i in other_points['y']:
        if y_i>=tip_displacement['y']:
            print('Y value out of bounds!')
    A = calculate_shape_coefficients_tracing(A0, other_points['y'], other_points['x'], N1, N2, chord = tip_displacement['y'], EndThickness = tip_displacement['x'])
    
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
    Au_P = [0.23993240191629417, 0.34468227138908186, 0.18125405377549103, 
            0.35371349126072665, 0.2440815012119143, 0.25724974995738387]
    Al_P = [0.18889012559339036, -0.24686758992053115, 0.077569769493868401,
            -0.547827192265256, -0.0047342206759065641, -0.23994805474814629]
    # NACA0012
    # Au_P =  [0.10887, 0.1187, 0.07843, 0.12084, 0.07919, 0.09840]
    # Al_P =  [0.11117, 0.1000, 0.1239, 0.06334, 0.11539, 0.10400] 
    # Passive shape coefficients for parent
    # Au_P = [.5,.4,.3]
    # Active shape coefficients for parent
    # Al_P = [.5,.1,.1]
   
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
    # AC_u1 = 0.1487            #Adimensional
    # AC_u2 = 0.10843          #Adimensional
    # AC_u3 = 0.15084                #Adimensional
    # AC_u4 = 0.10919              #Adimensional
    # AC_u5 = 0.12840             #Adimensional
    
    # Passive shape coefficients for child
    AC_u = [.25, .25, .25]
    
    # AC_u1 = 0.34468227138908186                #Adimensional
    # AC_u2 = 0.18125405377549103                 #Adimensional
    # AC_u3 = 0.35371349126072665                #Adimensional
    # AC_u4 = 0.2440815012119143                 #Adimensional
    # AC_u5 = 0.25724974995738387                 #Adimensional
    #Spar position for cruise (adiminesional because the chord will still be calculated)
    psi_spars =  [0.1, 0.2, 0.3]

    #==============================================================================
    # Calculate dependent coefficients
    #==============================================================================
    Au_C, Al_C, c_C, spar_thicknesses = calculate_dependent_shape_coefficients(
                                                        AC_u,
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
    plt.ylabel(r'$\zeta^p$', fontsize = 14)
    plt.ylim([-0.06,0.17])
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend(loc=1)
    plt.show()
    
    if morphing_direction == 'forwards':
        print('chords', c_P, c_C)
        # Calculate initial lengths
        strains, av_strains = calculate_strains(Au_P, Al_P, c_P, Au_C, Al_C, c_C, deltaz, psi_spars)
        
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
            print((length_children-length_parent)/length_parent)