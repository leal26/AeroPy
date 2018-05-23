# -*- coding: utf-8 -*-
"""
Objective: create an airfoil with a leading edge restriction, same upper length
restriction, othogonal upper spars and constant thicknesses in four places

Created on Mon Oct 17 10:36:34 2016

@author: Pedro
"""
import os
import math
import numpy as np
from numpy.linalg import inv

from aeropy.airfoil_module import CST
from aeropy.CST.module_3D import CST_3D
from aeropy.CST.module_2D import calculate_c_baseline, calculate_psi_goal, calculate_spar_direction

# Just as quick trick, to make upper morph I just mirror the image in regards to x
inverted = False
# Defines if basckwards or forwards morphing
morphing_direction = 'forwards'
    
#==============================================================================
# Calculate dependent shape function parameters
#==============================================================================
def calculate_dependent_shape_coefficients(BP_p, BA_p, BP_c, chord_p, sweep_p, 
                  twist_p, delta_TE_p, sweep_c, twist_c, eta_sampling, psi_spars, 
                  morphing='camber'):
    """Calculate  dependent shape coefficients for children configuration for a 4 order
    Bernstein polynomial and return the children upper, lower shape 
    coefficients, children chord and spar thicknesses. _P denotes parent parameters"""
    def calculate_BP_c0(BP_c0, c_P, deltaz, j):
        Au_C = extract_A(BP_c,j)
        Au_P = extract_A(BP_p,j)
        c_C = calculate_c_baseline(c_P, Au_C, Au_P, deltaz)
        BP_c[0][j] = np.sqrt(c_P/c_C)*Au_P[0]
        return BP_c[0][j], c_C
    
    def extract_A(B,j):
        # Extracting shape coefficient data for column j
        A = []
        for i in range(n+1):
            A.append(B[i][j])
        return A
    # Bersntein Polynomial
    def K(r,n):
        K=math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
        return K

    # Bernstein Polynomial orders (n is for psi, and m for eta)
    n = len(BP_p) - 1
    m = len(BP_p[0]) - 1
    p = len(psi_spars)
    q = len(eta_sampling)
    # print p,q,n,m    
    # Define chord, sweep, and twist functions for parent
    chord_p = CST(eta_sampling, chord_p['eta'][1], chord_p['initial'], Au=chord_p['A'], 
                  N1=chord_p['N1'], N2=chord_p['N2'], deltasLE=chord_p['final'])
    sweep_p = CST(eta_sampling, sweep_p['eta'][1], deltasz = sweep_p['final'], 
                  Au=sweep_p['A'], N1=sweep_p['N1'], N2=sweep_p['N2'])
    chord_p = chord_p[::-1]
    sweep_p = sweep_p
    twist_p = CST(eta_sampling, twist_p['eta'][1], twist_p['initial'], Au=twist_p['A'], 
                  N1=twist_p['N1'], N2=twist_p['N2'], deltasLE=twist_p['final'])
    delta_TE_p = CST(eta_sampling, delta_TE_p['eta'][1], delta_TE_p['initial'], 
                    Au=delta_TE_p['A'], N1=delta_TE_p['N1'], N2=delta_TE_p['N2'], 
                    deltasLE=delta_TE_p['final'])
    # Initialize chord, sweep, and twist functions for child             
    chord_c = []
    # Initialize child active matrix
    BA_c = []
    for i in range(n+1):
        temp = []
        for j in range(m+1):
            temp.append(0)
        BA_c.append(temp,)
    # Find upper shape coefficient though iterative method since Au_0 is unknown
    # via fixed point iteration
    for k in range(q):
        error = 9999
        BP_c0 = BP_p[0][k]
        while error > 1e-9:
            before = BP_c0
            c_P = chord_p[k]
            deltaz = delta_TE_p[k]
            [BP_c0,c_c] = calculate_BP_c0(BP_c0, c_P, deltaz, k)
            error = abs(BP_c0 - before)
        BP_c[0][k] = BP_c0
        BA_c[0][k] = np.sqrt(c_P/c_c)*BA_p[0][k]
        chord_c.append(c_c)
        print c_c, BP_c0, np.sqrt(c_P/c_c)*BA_p[0][k]
    # Calculate thicknessed and tensor C for the constraint linear system problem
    psi_A_c = []
    
    if morphing == 'camber':
        f = np.zeros((q,p))
        for l in range(q):
            # Converting everything from 3D to 2D framework
            Au_P = extract_A(BP_p,l)
            Al_P = extract_A(BA_p,l)
            Au_C = extract_A(BP_c,l)
            c_P = chord_p[l]
            c_C = chord_c[l]
            deltaz = delta_TE_p[l]
            
            # psi/xi coordinates for lower surface of the children configuration
            psi_lower_children = []
            xi_upper_children = []

            # psi_baseline, Au_baseline, Au_goal, deltaz, c_baseline, c_goal
            psi_upper_children = []
            for j in range(len(psi_spars)):
                psi_upper_children.append(calculate_psi_goal(psi_spars[j], Au_P, Au_C, deltaz,
                                       c_P, c_C))
            # Calculate xi for upper children. Do not care about lower so just gave it random shape coefficients
            xi_upper_children = CST(psi_upper_children, 1., deltasz= [deltaz/2./c_C, deltaz/2./c_C],  Al= Au_C, Au =Au_C)
            xi_upper_children = xi_upper_children['u']

            # print xi_upper_children
            
            #Debugging section
            x = np.linspace(0,1)
            y = CST(x, 1., deltasz= [deltaz/2./c_C, deltaz/2./c_C],  Al= Au_C, Au =Au_C)
            # plt.plot(x,y['u'])
            # plt.scatter(psi_upper_children, xi_upper_children)
            # plt.grid()
            # plt.show()
            # BREAK
            for k in range(len(psi_spars)):
                xi_parent = CST(psi_spars, 1., deltasz= [deltaz/2./c_P, deltaz/2./c_P],  Al= Al_P, Au =Au_P)
                delta_k_P = xi_parent['u'][k]-xi_parent['l'][k]
                t_k = c_P*(delta_k_P)
                # Claculate orientation for children
                s_k = calculate_spar_direction(psi_spars[k], Au_P, Au_C, deltaz, c_C)
                psi_l_k = psi_upper_children[k]-delta_k_P/c_C*s_k[0]
                xi_l_k = xi_upper_children[k]-delta_k_P/c_C*s_k[1]

                psi_lower_children.append(psi_l_k)

                f_y = 0
                for j in range(m+1):
                    f_y += BA_c[0][j]*(1-psi_l_k)**n*(K(j,m)*eta_sampling[l]**j*(1-eta_sampling[l])**(m-j))
                    # print j, eta_sampling[l]**j*(1-eta_sampling[l])**(m-j), BA_c[0][j]
                f[l][k] = (2*xi_l_k + psi_l_k*deltaz/c_C)/(2*(psi_l_k**0.5)*(psi_l_k-1))  - f_y
                
                # print 'f_y',f_y, eta_sampling[l],m,j
            # Store new children psi values
            psi_A_c.append(psi_lower_children)

        # Initialize F (avoiding using numpy)
        F = np.zeros([q,p,m+1,n])
        # F = []
        # for l in range(q):
            # tempk = []
            # for k in range(p):
                # tempj = []
                # for j in range(m+1):
                    # tempi = []
                    # for i in range(n):
                        # tempi.append(0.0)
                    # tempj.append(tempi)
                # tempk.append(tempj)
            # F.append(tempk)

        #j is the row dimension and i the column dimension in this case
        for l in range(q):
            for k in range(p):
                for j in range(m+1):
                    for i in range(n):
                        #Because in Python counting starts at 0, need to add 1 to be
                        #coherent for equations
                        ii = i + 1
                        Sx = K(ii,n)*(psi_A_c[l][k]**ii)*(1-psi_A_c[l][k])**(n-ii)
                        Sy = K(j,m)*(eta_sampling[l]**j)*(1-eta_sampling[l])**(m-j)
                        F[l][k][j][i] = Sx*Sy

        # print len(F), len(F[0]), len(F[0][0]), len(F[0][0][0])

        # Unfolding tensor
        F_matrix = np.zeros((n**2,n**2))
        f_vector = np.zeros((n**2,1))
        for l in range(q):
            for k in range(p):
                for j in range(m+1):
                    for i in range(n):
                        ii = n*(l) + k 
                        jj = n*(j) + i

                        F_matrix[ii][jj] = F[l][k][j][i]
                        f_vector[ii] = f[l][k]

        solution = np.linalg.solve(F_matrix,f_vector)

        for j in range(m+1):
            for i in range(n):
                jj = n*(j) + i 
                BA_c[i+1][j] = solution[jj][0]
        print BA_c
    return BA_c, chord_c
        
if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm

    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Inputs
    # One of the diameters
    initial_chord = 1.
    final_chord = 1.
    # Nosecone height
    span = 2.
    # Passive shape coefficients for parent
    P_p = [.5,.4,.3]
    # Active shape coefficients for parent
    A_p = [.5,.1,.1]
    # Passive shape coefficients for child
    P_c = [1.,.7]
    # location of the nosecone tip
    initial_nosecone_x = 0
    final_nosecone_x = 0.
    # Class coefficient for chord distribution (Nb=.5, elliptical, Nb=1, Haack series)
    Nb = 0.0
    # Axis of rotation
    axis = [final_nosecone_x - initial_nosecone_x, span,0]
    # Location of spars
    psi_spars = [0.3,0.7]
    # Sampling of spanwise dimension
    eta_sampling = [.0,1.]
    # for i in range(len(psi_spars)):
        # eta_sampling.append(i/(len(psi_spars)-1))
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BP_p = [[P_p[0],  P_p[0]],
            [P_p[1],  P_p[1]],
            [P_p[2],  P_p[2]]]
    BA_p = [[A_p[0],  A_p[0]],
            [A_p[1],  A_p[1]],
            [A_p[2],  A_p[2]]]
    BP_c = [[P_p[0],  P_p[0]], #Initial guess for leading edge is the parent
            [P_c[0],  P_c[0]],
            [P_c[1],  P_c[1]]]
    Na = 0.0
    mesh =(40,40)
    N={'eta':[0,1], 'N1':[.5,.5], 'N2':[1., 1.]}
    chord = {'eta':[0,1], 'A':[0.], 'N1':Na, 'N2':Nb, 
             'initial':initial_chord,'final':final_chord}
    sweep = {'eta':[0,1], 'A':[0.], 'N1':Nb, 'N2':Na, 
             'initial':initial_nosecone_x, 'final':final_nosecone_x}
    twist = {'eta':[0,1], 'A':[0.], 'N1':1, 'N2':1,'initial':0., 'final':0,
             'psi_root':0.,'zeta_root':0, 'axis':axis}
    delta_TE = {'eta':[0,1], 'A':[0.], 'N1':Nb, 'N2':Na, 
             'initial':0, 'final':0}
        
    BA_c, chord_c = calculate_dependent_shape_coefficients(
                     BP_p, BA_p, BP_c, chord, sweep, twist, delta_TE, sweep, twist, eta_sampling,
                     psi_spars, morphing='camber')
    chord_C = {'eta':[0,1], 'A':[0.], 'N1':Na, 'N2':Nb, 
             'initial':chord_c[0],'final':chord_c[1]}     

    Data = CST_3D(BP_c, BA_c,  span, N, mesh, chord_C, sweep, twist)
    X_u = Data['upper']['z']
    X_l = Data['lower']['z']
    Y_u = Data['upper']['y']
    Y_l = Data['lower']['y']
    Z_u = Data['upper']['z']
    Z_l = Data['lower']['z']
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # surf_u = ax.plot_surface(X_u, Z_u, Y_u, cmap=plt.get_cmap('jet'),
                       # linewidth=0, antialiased=False)
    # surf_l = ax.plot_surface(X_l, Z_l, Y_l, cmap=plt.get_cmap('jet'),
                       # linewidth=0, antialiased=False)
    
    # Customize the z axis.
    # ax.set_zlim(0, 4)

    max_range = np.array([X_u.max()-X_u.min(), X_l.max()-X_l.min(), Z_u.max()-Z_l.min(),
                          Y_u.max()-Y_u.min(), Y_l.max()-Y_l.min()]).max() / 2.0

    mid_x = (X_u.max()+X_l.min()) * 0.5
    mid_y = (Y_u.max()+Y_l.min()) * 0.5
    mid_z = (Z_u.max()+Z_l.min()) * 0.5
    # ax.set_xlim(mid_x - max_range, mid_x + max_range)
    # ax.set_ylim(mid_z - max_range, mid_z + max_range)
    # ax.set_zlim(mid_y - max_range, mid_y + max_range)
    # plt.xlabel('x')
    # plt.ylabel('z')
    # plt.show()

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot_trisurf(X_u.flatten(),  Y_u.flatten(), Z_u.flatten(), linewidth=0.2, antialiased=True,color='r', alpha=1.)
    ax.plot_trisurf(X_l.flatten(),  Y_l.flatten(), Z_l.flatten(), linewidth=0.2, antialiased=True,color='r', alpha=1.)
    
    # Doing for the original configuration
    Data = CST_3D(BP_p, BA_p,  span, N, mesh, chord_C, sweep, twist)
    X_u = Data['upper']['z']
    X_l = Data['lower']['z']
    Y_u = Data['upper']['y']
    Y_l = Data['lower']['y']
    Z_u = Data['upper']['z']
    Z_l = Data['lower']['z']
    ax.plot_trisurf(X_u.flatten(),  Y_u.flatten(), Z_u.flatten(), linewidth=0.2, antialiased=True,color='b', alpha=1.)
    ax.plot_trisurf(X_l.flatten(),  Y_l.flatten(), Z_l.flatten(), linewidth=0.2, antialiased=True,color='b', alpha=1.)
    # ax.plot([sweep['initial'] + twist['psi_root']*chord['initial'],sweep['initial'] + twist['psi_root']*chord['initial'] + axis[0]],[0,axis[1]])
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    BREAK
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
            s = calculate_spar_direction(psi_i, Au_C, Au_P, deltaz, c_P, spar_thicknesses)

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