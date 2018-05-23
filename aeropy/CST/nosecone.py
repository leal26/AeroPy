from __future__ import print_function
from scipy.interpolate import interp1d
import numpy as np
import math

from aeropy.airfoil_module import CST

def taper_function(eta, shape = 'linear', points = {'eta':[0,1], 'chord':[1,.7]}):
    """Calculate chord along span of the wing.

    - If linear, taper function is a conjuction of lines connecting points
    - Possible shapes are the same as interp1d: ('linear', 'nearest',
        'zero', 'slinear', 'quadratic', 'cubic' where 'zero', 'slinear',
         'quadratic' and 'cubic' refer to a spline interpolation of zeroth,
          first, second or third order"""
    function = interp1d(points['eta'], points['chord'])
    return function(eta)

def twist_function(eta, shape = 'linear', points = {'eta':[0,1], 'delta_twist':[0,.1]}):
    """Calculate chord along span of the wing.

    - If linear, taper function is a conjuction of lines connecting points
    - Possible shapes are the same as interp1d: ('linear', 'nearest',
        'zero', 'slinear', 'quadratic', 'cubic' where 'zero', 'slinear',
         'quadratic' and 'cubic' refer to a spline interpolation of zeroth,
          first, second or third order"""
    function = interp1d(points['eta'], points['delta_twist'])
    return function(eta)

def CST_3D(Bu, Bl, span, N={'eta':[0,1], 'N1':[.5, .5], 'N2':[1., 1.], 'chord':[1., 0]},
           mesh = (100,100), chord = {'eta':[0,1], 'A':[1.], 'N1':1, 'N2':1, 'initial_chord':1.}, 
           sweep = {'eta':[0,1], 'A':[1.], 'N1':1, 'N2':1, 'x_LE_initial':0, 'x_LE_final':0}):
    """
    - Bu: upper shape coefficients
    - Bl: lower shape coefficients
    - mesh: list of number of points in x and y
    """
    def S(B, psi, eta):
        """ Cross section shape function. Validated for high dimensions.
           To debug just verify if it turns all ones when B=ones"""
        def S_i(r, n, psi):
            """Shape function"""
            value = K(r,n)*(psi**r)*(1.-psi)**(n-r)
            return value

        # Bersntein Polynomial
        def K(r,n):
            K=math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
            return K

        Nx = len(B)-1
        Ny = len(B[0])-1

        output = 0
        for i in range(Nx+1):
            for j in range(Ny+1):
                output += B[i][j]*S_i(i, Nx, psi)*S_i(j, Ny, eta)
        return output

    def C(N, psi, eta):
        """Class function"""
        N1 = interp1d(N['eta'], N['N1'])
        N2 = interp1d(N['eta'], N['N2'])
        output = ((psi)**N1(eta))*((1.-psi)**N2(eta))
        return output

    psi = np.linspace(0,1,mesh[0])
    eta = np.linspace(0,1,mesh[1])

    zeta_u = np.zeros(mesh)
    zeta_l = np.zeros(mesh)
    for i in range(mesh[0]):
        for j in range(mesh[1]):
            zeta_u[j][i] = C(N, psi[i], eta[j])*S(Bu, psi[i], eta[j])
            zeta_l[j][i] = -C(N, psi[i], eta[j])*S(Bl, psi[i], eta[j])
    print(eta)
    print(chord['initial_chord'])
    print(chord['A'])
    print(chord['N1'], chord['N2'])
    chord_distribution = CST(eta, chord['eta'][1], chord['initial_chord'], Au=chord['A'], N1=chord['N1'], N2=chord['N2'])
    sweep_distribution = CST(eta, sweep['eta'][1], deltasz = sweep['x_LE_final']-.5*chord['initial_chord'], Au=sweep['A'], N1=sweep['N1'], N2=sweep['N2'])
    chord_distribution = chord_distribution[::-1]
    sweep_distribution = sweep_distribution
    # taper_function(eta, shape = 'linear', N)
    x = np.zeros(len(psi))
    for i in range(len(x)):
        x[i] = psi[i]*chord_distribution[i]
    print(chord_distribution)
    print(sweep_distribution)
    print(x)
    print(psi)
    y = eta

    X = np.zeros(mesh)
    Y = np.zeros(mesh)
    Z_u = np.zeros(mesh)
    Z_l = np.zeros(mesh)
    for i in range(mesh[0]):
        for j in range(mesh[1]):
            X[j][i] = psi[i]*chord_distribution[j] - sweep_distribution[j] -.5*chord['initial_chord']
        
            Y[j][i] = span*eta[j]
            Z_u[j][i] = zeta_u[j][i]*chord_distribution[j]
            Z_l[j][i] = zeta_l[j][i]*chord_distribution[j]
    return [X,Y,Z_u,Z_l]

if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm

    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Inputs
    # One of the diameters
    initial_chord = 1.
    # Nosecone height
    span = 4.
    # Shape coefficient for cross section (if A=1, circular, otherwise it is an ellipse)
    A = 1.
    # location of the nosecone tip
    nosecone_x = 0.2
    # Class coefficient for chord distribution (Nb=.5, elliptical, Nb=1, Haack series)
    Nb = 1.
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #B = [[1,1], [1.,1]]
    B = [[A], [A]]
    Na = 1.
    x = np.linspace(0,1)
    
    [X,Y,Z_u, Z_l] = CST_3D(B, B, mesh =(50,50), span=span,
                            N={'eta':[0,1], 'N1':[.5, .5], 'N2':[.5, .5]},
                            chord = {'eta':[0,1], 'A':[1.], 'N1':Na, 'N2':Nb, 'initial_chord':initial_chord},
                            sweep = {'eta':[0,1], 'A':[.5], 'N1':Nb, 'N2':Na, 'x_LE_final':nosecone_x})

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf_u = ax.plot_surface(X, Z_u, Y, cmap=plt.get_cmap('jet'),
                       linewidth=0, antialiased=False)
    surf_l = ax.plot_surface(X, Z_l, Y, cmap=plt.get_cmap('jet'),
                       linewidth=0, antialiased=False)
    # cset = ax.contour(X, Z_u, Y, zdir='z', offset=0, cmap=cm.coolwarm)
    # cset = ax.contour(X, Z_l, Y, zdir='z', offset=0,  cmap=cm.coolwarm)
    # cset = ax.contour(X, Z_u, Y, zdir='x', offset=-.1, cmap=cm.coolwarm)
    # cset = ax.contour(X, Z_l, Y, zdir='x', offset=-.1, cmap=cm.coolwarm)
    # cset = ax.contour(X, Z_u, Y, zdir='y', offset =0.5,  cmap=cm.coolwarm)
    # cset = ax.contour(X, Z_l, Y, zdir='y', offset =0.5,  cmap=cm.coolwarm)
    
    # Customize the z axis.
    ax.set_zlim(0, 4)

    max_range = np.array([X.max()-X.min(),  Z_u.max()-Z_l.min(), Y.max()-Y.min()]).max() / 2.0

    mid_x = (X.max()+X.min()) * 0.5
    mid_y = (Y.max()+Y.min()) * 0.5
    mid_z = (Z_u.max()+Z_l.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_z - max_range, mid_z + max_range)
    ax.set_zlim(mid_y - max_range, mid_y + max_range)
    plt.show()
