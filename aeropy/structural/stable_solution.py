import numpy as np
from scipy.optimize import minimize


class mesh_1D():
    def __init__(self, alpha=[1, 1], alpha_nodes=[0, 1], mesh_n=10, x2=0):
        self.n = mesh_n
        self.x_p = np.linspace(0, 1, mesh_n)
        self.dx_p = 1./(mesh_n-1)
        self.alpha = alpha
        self.alpha_nodes = alpha_nodes
        self.x2 = x2

        if self.alpha_nodes[0] != 0 or self.alpha_nodes[-1] != 1:
            raise Exception('Need to define alpha_x for whole domain')

        if len(self.alpha_nodes) != len(alpha):
            raise Exception('Number of alpha nodes and values must be same')
        if list(self.alpha_nodes) != sorted(self.alpha_nodes):
            raise Exception('Alpha nodes must be in increasing order')

        self.nodes_xc = self.mesh_child()

    def mesh_child(self):
        self.n = len(self.x_p)
        x = self.x_p
        counter = 0
        self.x_c = self.x_p.copy()
        self.alpha_x = []
        for i in range(len(self.alpha_nodes)-1):

            alpha_i = self.alpha_nodes[i]
            alpha_f = self.alpha_nodes[i+1]
            x_filtered = x[x >= alpha_i]
            x_filtered = x_filtered[x_filtered < alpha_f]
            self.x_c[counter:counter+len(x_filtered)] = self.alpha[i]*x_filtered
            counter += len(x_filtered)
            self.alpha_x += len(x_filtered)*[self.alpha[i]]
        # Add last
        self.alpha_x.append(self.alpha[-1])
        self.alpha_x = np.array(self.alpha_x)
        self.x_c[-1] *= self.alpha[-1]


class properties():
    def __init__(self, young=70e9, poisson=.3, dimensions=[0.01, 0.01],
                 crosssection='square'):
        """Stores calculates material or structural properties. If 'square',
        dimensions = [width, height]"""
        self.young = young
        self.poisson = poisson
        self.dimensions = dimensions
        self.length = 1

        if crosssection == 'square':
            self.area = dimensions[0]*dimensions[1]
            self.inertia = (dimensions[0]*dimensions[1]**3)/12.
        if young < 0 or poisson < 0 or poisson > 1 or self.area <= 0:
            raise Exception('Material properties need to make sense')

class loads():
    def __init__(self, concentrated_load=np.array([[0, -1], ]), load_s=[1],
                       distributed_function = 0):
        self.concentrated_load = concentrated_load
        self.concentrated_s = load_s
        self.distributed_function = distributed_function


class boundary_conditions():
    def __init__(self, concentrated_load=np.array([[0, 0, 0], ]), load_x=[1],
                       distributed_load = 0):
        self.concentrated_load = concentrated_load
        self.concentrated_x = load_x
        self.concentrated_n = len(concentrated_load)
        self.distributed_load = None

        self.distributed_load = distributed_load

        if len(concentrated_load) != len(load_x):
            raise Exception('load values and x lists have to match')

class euler_bernoulle():
    def __init__(self, properties, load, load_type, geometry,
                 x = np.linspace(0,1), boundaries = {'x':[], 'u':[]}):
        self.properties = properties
        self.load = load
        self.load_type = load_type
        self.g = geometry
        self.g.x1_grid = x
        self.boundaries = boundaries

    def analytical_solutions(self):
        if self.load_type == 'concentrated':
            bp = self.properties
            self.g.D = self.load/(6*bp.young*bp.inertia) * \
                            np.array([0, 0, 3, -1, 0])
        elif self.load_type == 'distributed':
            bp = self.properties
            c2 = 6*(bp.length**2)*self.load/(24*bp.young*bp.inertia)
            c3 = -4*(bp.length)*self.load/(24*bp.young*bp.inertia)
            c4 = (1)*self.load/(24*bp.young*bp.inertia)
            self.g.D = np.array([0, 0, c2, c3, c4])

    def bending_strain(self):
        print(self.g.x1_grid)
        self.B = self.g.x3(self.g.x1_grid, 'x11')

    def free_energy(self):
        bp = self.properties
        ddu = self.g.x3(self.g.x1_grid, 'x11')
        self.phi = bp.young*bp.inertia/2*ddu**2

    def strain_energy(self):
        self.U = np.trapz(self.phi, self.g.x1_grid)

    def work(self):
        u = self.g.x3(self.g.x1_grid)
        if self.load_type == 'concentrated':
            self.W = self.load*u[-1]
        elif self.load_type == 'distributed':
            self.W = np.trapz(self.load*u, self.g.x1_grid)

    def residual(self):
        self.R = self.U - self.W

    def minimum_potential(self, x0=[0,0,0]):
        def to_optimize(x):
            self.g.D = [0,0] + list(x)
            self.work()
            self.free_energy()
            self.strain_energy()
            self.residual()
            return self.R

        # With bounds
        bounds = np.array(((-0.01,0.01),)*len(x0))

        res = minimize(to_optimize, x0, bounds=bounds)
        self.g.D = [0,0] + list(res.x)
        self.work()
        self.free_energy()
        self.strain_energy()
        self.residual()
        return(res.x, res.fun)

class structure():
    def __init__(self, geometry_parent, geometry_child, mesh, properties,
                 bc, model='beam'):
        # Defining geometries
        self.g_p = geometry_parent
        self.g_c = geometry_child

        self.model = model
        self.bc = bc
        self.mesh = mesh
        self.properties = properties

        # storage varaibles
        self.opt_x = []
        self.opt_f = []

    def u(self, input=None, diff=None):
        # If for a one time run, run for new input and revert back to original
        # input
        if input is not None:
            stored_x_p = self.mesh.x_p
            self.mesh.x_p = input
            self.mesh.mesh_child()

        parent = self.g_p.r(input=self.mesh.x_p, x2=self.mesh.x2,
                            input_type='x1', diff=diff)
        child = self.g_c.r(input=self.mesh.x_c, x2=self.mesh.x2,
                           input_type='x1', diff=diff)

        self.cosine_direction(diff=None)
        # Taking into consideration extension of the beam
        child[0] *= self.mesh.alpha_x
        output = child - parent

        for i in range(self.mesh.n):
            output[:, i] = np.matmul(self.R[i], output[:, i])
        if diff is 'x1':
            dR = self.cosine_direction(diff='x1')
            parent = self.g_p.r(input=self.mesh.x_p, x2=self.mesh.x2,
                                input_type='x1', diff=None)
            child = self.g_c.r(input=self.mesh.x_c, x2=self.mesh.x2,
                               input_type='x1', diff=None)
            position_delta = child - parent
            output[:, i] += np.matmul(dR[i], position_delta[:, i])

        if input is not None:
            self.mesh.x_p = stored_x_p
            self.mesh.mesh_child()
        return(output)

    def calculate_position(self, input=None, diff=None):
        # If for a one time run, run for new input and revert back to original
        # input
        if input is not None:
            stored_x_p = self.mesh.x_p
            self.mesh.x_p = input
            self.mesh.mesh_child()

        self.r_p = self.g_p.r(input=self.mesh.x_p, x2=self.mesh.x2,
                              input_type='x1', diff=diff)
        self.r_c = self.g_c.r(input=self.mesh.x_c, x2=self.mesh.x2,
                              input_type='x1', diff=diff)
        if diff is not None:
            r_p, r_c = self.r_p, self.r_c
            self.calculate_position(input=input, diff=None)
            return(r_p, r_c)
        if input is not None:
            self.mesh.x_p = stored_x_p
            self.mesh.mesh_child()

    def cosine_direction(self, diff=None, g=None):
        """Calculate cosine matrix between a rectangular cartesian system and
        a curvilinear(or given) coordinate system. Returns matrix with shape
        (number of nodes, Ri, Rj). Different from rest of code where (x,y,n)"""
        def dot(a, b):
            a_1, a_2 = a
            b_1, b_2 = b
            return(a_1*b_1 + a_2*b_2)
        self.R = np.zeros((self.mesh.n, 2, 2))
        e = np.eye(2)
        for k in range(self.mesh.n):
            for i in range(2):
                for j in range(2):
                    gi = e[i]
                    if g is None:
                        gj = self.g_p.g(j+1, np.array([self.mesh.x_c[k]]),
                                        diff=diff)
                    else:
                        gj = g[j]

                    self.R[k][i][j] = dot(gi, gj)
        return(self.R)

    def uij(self, i, j, diff=None, input_type='x1'):
        '''Indexes here are from 1 to n. So +=1 compared to rest'''
        # TODO: makes this more optimal(calculating u multiple times)

        ui_j = self.u(diff='x%i' % (j))[i-1]
        ui = self.u()[i-1]

        for l in range(1, 3):
            ui_j += self.g_c.christoffel(i, j, l, self.mesh.x_p,
                                         self.mesh.x2)*ui

        return(ui_j)

    def strain(self):
        # For cartesian
        self.epsilon = np.zeros([2, 2, self.mesh.n])
        for i in range(2):
            for j in range(2):
                ii = i + 1
                jj = j + 1
                self.epsilon[i][j] = .5*(self.uij(ii, jj) +
                                         self.uij(jj, ii))
                # for m in range(2):
                #     mm = m + 1
                #     self.epsilon[i][j] -= .5*(self.uij(mm, jj) *
                #                              self.uij(mm, ii))
        return(self.epsilon)

    def stress(self, loading_condition='uniaxial'):
        E = self.properties.young
        nu = self.properties.poisson
        if loading_condition == 'uniaxial':
            self.sigma = E*self.epsilon
        if loading_condition == 'plane_stress':
            self.sigma = E/(1-nu**2)*(1-nu)*self.epsilon
            for k in range(2):
                self.sigma[k][k] += E/(1-nu**2)*nu*self.epsilon[k][k]
        return(self.sigma)

    def strain_energy(self):
        energy = 0
        for i in range(len(self.sigma)):
            for j in range(len(self.sigma[i])):
                for k in range(len(self.sigma[i][j])):
                    if k == 0 or k == self.mesh.n - 1:
                        multiplier = .5*self.properties.area*self.mesh.dx_p/2.
                    else:
                        multiplier = .5*self.properties.area*self.mesh.dx_p
                    energy += multiplier*self.sigma[i][j][k]*self.epsilon[i][j][k]
        return(energy)

    def work(self):
        energy = 0
        for i in range(self.bc.concentrated_n):
            u = self.u(np.array([self.bc.concentrated_x[i]]))
            for j in range(2):
                energy += self.bc.concentrated_load[i][j] * u[j][0]
        return(energy)

    def residual(self, input=None, input_type = 'Strain',
                loading_condition = None, input_function = None):
        if input is not None:
            if input_function is not None:
                input = input_function(input)
            if input_type is 'Strain':
                self.mesh.alpha = 1+np.array(input)
                self.mesh.mesh_child()
            elif input_type is 'Geometry':
                self.g_c.a = input
                self.mesh.mesh_child()
                self.calculate_position()
            self.strain()
            self.stress(loading_condition = loading_condition)
        energy = self.strain_energy() - self.work()
        return(energy)

    def find_stable(self, x0=[0], bounds=None, input_type = 'Strain',
                    loading_condition = 'uniaxial', input_function = lambda x:x):
        def _callback(x):
            input = (bounds[:,1]-bounds[:,0])*x + bounds[:,0]
            self.opt_x.append(input)
            self.opt_f.append(self.residual(input, input_type = input_type,
                                 loading_condition = loading_condition,
                                 input_function = input_function))
        def to_optimize(x):
            input = (bounds[:,1]-bounds[:,0])*x + bounds[:,0]
            output = self.residual(input, input_type = input_type,
                                 loading_condition = loading_condition,
                                 input_function = input_function)

            # print(output, input)
            return output

        # With bounds
        try:
            x0_nd = (x0-bounds[:,0]) / (bounds[:,1]-bounds[:,0])
            bounds_nd = np.array([[-1,1],]*len(x0))

            self.opt_x = [x0]
            self.opt_f = [to_optimize(x0_nd)]

            res = minimize(to_optimize, x0_nd, bounds=bounds_nd, callback=_callback)
        except(TypeError):
            bounds = np.array(((0,1),)*len(x0))
            self.opt_x = [x0]
            self.opt_f = [to_optimize(x0)]

            res = minimize(to_optimize, x0, callback=_callback)
        x = (bounds[:,1]-bounds[:,0])*res.x + bounds[:,0]
        return(x, res.fun)

    def sweep_strains(self, strains, strains_x, reorder=None,
                      loading_condition = 'uniaxial'):
        energy_list = []
        residual_list = []
        n = len(strains)
        for i in range(n):
            self.mesh.alpha = 1 + strains[i]
            self.mesh.alpha_x = strains_x
            self.mesh.mesh_child()
            self.strain()
            self.stress()
            residual_list.append(self.residual(input_type = 'Strain',
                        loading_condition = loading_condition))
            energy_list.append(self.strain_energy())
        if reorder is not None:
            residual_list = np.resize(residual_list, reorder)
            energy_list = np.resize(energy_list, reorder)
        return(energy_list, residual_list)

    def sweep_geometries(self, geom_variables, input_function, reorder=None,
                         loading_condition = 'plane_stress'):
        energy_list = []
        residual_list = []
        n = len(geom_variables)
        for i in range(n):
            print(i)
            input = geom_variables[i]
            residual_list.append(self.residual(input, input_type = 'Geometry',
                                               input_function = input_function,
                                               loading_condition = loading_condition))
            energy_list.append(self.strain_energy())
        if reorder is not None:
            residual_list = np.resize(residual_list, reorder)
            energy_list = np.resize(energy_list, reorder)
        return(energy_list, residual_list)
