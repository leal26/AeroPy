import numpy as np


class structure():
    def __init__(self, geometry_parent, geometry_child, model='beam',
                 young=70e9, poisson=.3, area=0.001):
        self.g_p = geometry_parent
        self.g_c = geometry_child
        self.model = model
        self.a_p = [0, 0, 0, 0]
        self.young = young
        self.poisson = poisson
        self.area = area

    def u(self, input, x2=0, diff=None, input_type='x1'):
        parent = self.g_p.r(input, x2=x2, input_type=input_type, diff=diff)
        child = self.g_c.r(input, x2=x2, input_type=input_type, diff=diff)
        output = child - parent
        return(output)

    def uij(self, i, j, input, x2=0, diff=None, input_type='x1'):
        '''Indexes here are from 1 to n. So +=1 compared to rest'''
        ui_j = self.u(input, diff='x%i' % (j))[i-1]
        ui = self.u(input)[i-1]
        for l in range(1, 3):
            ui_j += self.g_c.christoffel(i, j, l, input, x2)*ui
        return(ui_j)

    def strain(self, input, x2=0):
        # For cartesian
        self.epsilon = np.zeros([2, 2, len(input)])
        for i in range(2):
            for j in range(2):
                ii = i + 1
                jj = j + 1
                self.epsilon[i][j] = .5*(self.uij(ii, jj, input) +
                                         self.uij(jj, ii, input))
        # christoffel_122 = self.g_c.christoffel(input, x2)
        return(self.epsilon)

    def stress(self):
        self.sigma = self.young/(1+self.poisson)*self.epsilon
        # for main diagonal components
        const = self.young/(1+self.poisson) * self.poisson/(1 - 2*self.poisson)
        for i in range(2):
            for k in range(2):
                self.sigma[i][i] += const*self.epsilon[k][k]
        return(self.sigma)

    def strain_energy(self, z1):
        self.strain(z1)
        self.stress()
        return(np.tensordeot(self.sigma, self.epsilon))
