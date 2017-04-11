from scipy.interpolate import interp1d
import numpy as np
import math

def taper_function(eta, shape = 'linear', points = {'eta':[0,1], 'c':[1,.7]}):
    """Calculate chord along span of the wing.

    - If linear, taper function is a conjuction of lines connecting points
    - Possible shapes are the same as interp1d: ('linear', 'nearest',
        'zero', 'slinear', 'quadratic', 'cubic' where 'zero', 'slinear',
         'quadratic' and 'cubic' refer to a spline interpolation of zeroth,
          first, second or third order"""
    function = interp1d(points['eta'], points['c'])
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

def CST_3D(Bu, Bl, N={'eta':[0,1], 'N1':[.5, .5], 'N2':[1., 1.]},
           mesh = (100,100)):
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
            zeta_u[i][j] = C(N, psi[i], eta[j])*S(Bu, psi[i], eta[j])
            zeta_l[i][j] = -C(N, psi[i], eta[j])*S(Bl, psi[i], eta[j])
    x = psi
    y = eta
    z_u = zeta_u
    z_l = zeta_l
    return [x,y,z_u,z_l]

if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm

    B = [[1,1], [1.,1]]

    [x,y,Z_u, Z_l] = CST_3D(B, B, mesh =(50,50),
                            N={'eta':[0,1], 'N1':[.5, .5], 'N2':[1., 1.]})
    X,Y = np.meshgrid(x,y)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf_u = ax.plot_surface(X, Y, Z_u, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    surf_l = ax.plot_surface(X, Y, Z_l, cmap=cm.coolwarm_r,
                       linewidth=0, antialiased=False)
    # Customize the z axis.
    ax.set_zlim(-0.5, .5)
    plt.show()
