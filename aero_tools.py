# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 17:26:19 2015

@author: Pedro Leal
"""
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
    
if __name__ == '__main__':
    print LLT_calculator(alpha_L_0_root=1., c_D_xfoil=0.01)
