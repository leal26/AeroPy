import math
import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize

from aeropy.geometry.airfoil import CST
from aeropy.CST_2D import (dxi_u, ddxi_u, calculate_c_baseline, K,
                           calculate_psi_goal, calculate_spar_direction, S)


class CoordinateSystem(object):
    """class for a polynomial function
    """

    def __init__(self, D, **kwargs):
        # Other inputs vary according to shape class
        for key, value in kwargs.items():
            setattr(self, key, value)
        # D is always the shape parameters
        self.D = D

    @property
    def D(self):
        return self._D

    @D.setter
    def D(self, values):
        if hasattr(self, 'name') and self.name == 'pCST':
            self._D = values
            # Update internal variables if class already populated
            if hasattr(self, 'cst'):
                bool, given_length, required_length = self._check_input(values)
                if bool:
                    self._pCST_update()
                else:
                    raise Exception('Input length is incorrect. ' +
                                    'It is %i and it should be %i' % (
                                        given_length, required_length))

        else:
            if self.n == 1 and (isinstance(values, float) or isinstance(values, np.float64)):
                self._D = values
            else:
                if len(values) != self.n:
                    self._D = np.zeros(self.n)
                    self._D[:len(values)] = values
                else:
                    self._D = values

    @property
    def x1_grid(self):
        return self._x1_grid

    @x1_grid.setter
    def x1_grid(self, values):
        if type(values) == list:
            values = np.array(values)

        if self.name == 'pCST':
            for i in range(self.p):
                if i != 0 or not self.rigid_LE:
                    self.cst[i].x1_grid = values[self.cst[i].indexes]
        self._x1_grid = values

    def _x1(self, x1, diff=None, offset=True):
        if offset:
            x1 = x1 - self.offset_x
        if diff is None:
            return x1
        elif diff == 'x1':
            try:
                return np.ones(len(x1))
            except:
                return 1
        elif diff == 'x11':
            try:
                return np.zeros(len(x1))
            except:
                return 0
        elif diff == 'theta3':
            return(self.a[2, :, 0])
        elif diff == 'theta1':
            output = 1/np.sqrt(1 + self.x3(x1, 'x1', False)**2)
            output[np.isnan(output)] = 0
            return output
            # return np.ones(len(x1))
        elif diff == 'theta11':
            # return -self.x1(x1, 'theta1')**4*self.x3(x1, 'x1')*self.x3(x1, 'x11')
            # dr = self.r(x1, diff='x1')
            # ddr = self.r(x1, diff='x11')
            # a1 = 1/np.sqrt(np.einsum('ij,ij->i',dr, dr))**3
            # a2 = np.einsum('ij,ij->i',dr, ddr)
            # return np.multiply(a1, a2)*self.x1(x1, 'theta1')
            # return np.zeros(len(x1))
            return -self.x1(x1, 'theta1', False)**4*self.x3(x1, 'x1', False)*self.x3(x1, 'x11', False)
        elif diff == 'theta31' or diff == 'theta13':
            return(-self.x3(x1, 'theta11', False))
            # return(np.zeros(len(x1)))
        else:
            return(np.zeros(len(x1)))

    def _x1_pCST(self, x1, diff=None):
        output = np.zeros(len(x1))
        if len(x1) == 1:
            c_min = 0
            for i in range(self.p):
                c_max = c_min + self.cst[i].chord
                if (x1[0] >= (c_min-self.tol)) & (x1[0] <= (c_max+self.tol)):
                    output[0] = self.cst[i].x1(x1, diff=diff)
                    break
                c_min = c_max
        else:
            for i in range(self.p):
                output[self.cst[i].indexes] = self.cst[i].x1(x1[self.cst[i].indexes], diff=diff)
        return output

    def x2(self, x1, x2_value=0, diff=None):
        if diff is None:
            return np.zeros(len(x1))
        if diff is 'x2' or diff is 'theta2':
            return np.ones(len(x1))

    @classmethod
    def polynomial(cls, D, chord=1, color='b', n=6, tol=1e-6):
        c = cls(D, chord=chord, color=color, n=n, tol=tol, offset_x=0)
        c.x1 = c._x1
        c.x3 = c._x3_poly
        c.name = 'polynomial'
        return c

    @classmethod
    def CST(cls, D, chord, color, N1=.5, N2=1.0, deltaz=0, tol=1e-6, offset=0,
            deltazLE=0, offset_x=0, origin=0):
        c = cls(D, chord=chord, color=color, n=len(D), N1=N1, N2=N2,
                deltaz=deltaz, tol=tol, offset=offset, deltazLE=deltazLE,
                offset_x=offset_x, origin=origin)
        c.x1 = c._x1
        c.x3 = c._x3_CST
        c.name = 'CST'
        c.zetaT = c.deltaz/c.chord
        c.zetaL = c.deltazLE/c.chord
        c.length = c.arclength()
        return c

    @classmethod
    def pCST(cls, D, chord=np.array([.5, .5]), color='b', N1=[1, 1.], N2=[1.0, 1.0],
             tol=1e-6, offset=0, offset_x=0, continuity='C2', free_end=False,
             root_fixed=False, dependent=None, length_preserving=False,
             rigid_LE=False):
        c = cls(D, chord=chord, color=color, n=len(D), N1=N1, N2=N2,
                tol=tol, offset=offset, offset_x=offset_x, continuity=continuity,
                free_end=free_end, root_fixed=root_fixed,
                length_preserving=length_preserving, rigid_LE=rigid_LE)
        c.x1 = c._x1_pCST
        c.x3 = c._x3_pCST
        c.name = 'pCST'
        c.p = len(chord)
        if dependent is None:
            c.dependent = len(N1)*[False]
        else:
            c.dependent = dependent

        # For conevinence n is order of CST and nn is number of shape coeff (n+1)
        if continuity == 'C2':
            c.nn = (len(D)-2)/c.p
        elif continuity == 'C1':
            c.nn = (len(D)-1)/c.p
        if (c.nn % 1) != 0:
            raise ValueError('Incorrect number of inputs D. Got nn=%f' % c.nn)
        else:
            c.nn = int(c.nn)
            if continuity == 'C1':
                c.n = c.nn - 1
            else:
                c.n = c.nn

        c.cst = []
        c.cst_p = []
        c.zetaL = []
        c.zetaT = []
        c.A0 = []
        c.A1 = c.p*[0]
        offset_s = 0
        offset_x = 0
        for i in range(c.p):
            j = i - 1
            # From shape coefficients 1 to n
            if continuity == 'C2':
                Ai = D[1+i*c.nn:1+(i+1)*c.nn]
            elif continuity == 'C1':
                Ai = D[i*c.nn:(i+1)*c.nn]
            if i == 0:
                c.A0.append(D[0])
                if c.root_fixed and c.N1[i] == 1:
                    c.zetaT.append(-D[0])
                else:
                    c.zetaT.append(D[-1])
                c.zetaL.append(0)
            else:
                if N1[i] == 1. and N2[i] == 1.:
                    offset_x += chord[j]
                    if continuity == 'C2':
                        ddj = c.n*c.cst[j].D[-3] - (N1[j]+c.n)*c.cst[j].D[-2]
                        c.A0.append((-chord[i]/chord[j]*ddj+Ai[0]*c.n)/(c.n+1))
                    elif continuity == 'C1':
                        c.A0.append(Ai[0])
                    c.zetaL.append(chord[j]/chord[i]*c.zetaT[j])
                    c.zetaT.append(-c.cst[j].D[-2] + c.zetaT[j] -
                                   c.A0[-1] + c.zetaL[i] - c.zetaL[j])
                else:
                    raise(NotImplementedError)

            if continuity == 'C2':
                Di = [c.A0[i]] + list(Ai) + [c.zetaT[i]]
            elif continuity == 'C1':
                Di = list(Ai) + [c.zetaT[i]]

            c.cst.append(CoordinateSystem.CST(Di, chord[i],
                                              color[i], N1=N1[i], N2=N2[i],
                                              deltaz=c.zetaT[-1]*chord[i],
                                              deltazLE=c.zetaL[-1]*chord[i],
                                              offset_x=offset_x, offset=offset))
            c.cst_p.append(copy.deepcopy(c.cst[-1]))
            c.cst[i].offset_s = offset_s
            offset_s += c.cst[i].length
        c.total_chord = sum([c.cst[i].chord for i in range(c.p)])
        c.total_length = sum([c.cst[i].length for i in range(c.p)])
        c.joint_s = sum([c.cst[i].offset_s for i in range(1, c.p)])
        return c

    def _pCST_update(self):
        offset_s = 0
        self.n_start = 0
        self.spar_i = 0

        for i in range(self.p):
            if i != 0 or not self.rigid_LE:
                if self.dependent[i]:
                    if not hasattr(self, 'g_independent'):
                        raise Exception('No geometry defined as reference for dependent')
                    self._calculate_spar(i)
                    self._update_dependent(i)
                else:
                    self._update_independent(i)
                self.cst[i].offset_s = offset_s
            offset_s += self.cst[i].length
        self.total_chord = sum([self.cst[i].chord for i in range(self.p)])

    def _update_independent(self, i):
        j = i - 1
        Ai0, Ai = self._select_D(i)
        # print('independent', Ai0, Ai)
        error = 999
        offset_x = 0
        while error > 1e-6:
            prev = np.array([self.cst[i].chord, self.cst[i].zetaT,
                             self.cst[i].zetaL, self.cst[i].D[0]])

            if i == 0:
                self.A0[i] = self.D[0]
                if self.root_fixed and self.N1[i] == 1:
                    self.zetaT[i] = -self.D[0]
                else:
                    self.zetaT[i] = self.D[-1]
                self.zetaL[i] = 0
            elif self.N1[i] == 1. and self.N2[i] == 1.:
                offset_x = self.cst[j].offset_x + self.cst[j].chord
                if self.continuity == 'C2':
                    ddj = self.n*self.cst[j].D[-3] - (self.N1[j]+self.n)*self.cst[j].D[-2]
                    self.A0[i] = (-self.cst[i].chord/self.cst[j].chord *
                                  ddj+Ai[0]*self.n)/(self.n+1)

                elif self.continuity == 'C1':
                    self.A0[i] = Ai[0]

                self.zetaL[i] = self.cst[j].chord/self.cst[i].chord*self.zetaT[j]
                self.zetaT[i] = -self.cst[j].D[-2] + self.zetaT[j] - self.A0[i] + \
                    self.zetaL[i] - self.zetaL[j]
            else:
                raise(NotImplementedError)

            An = self._calculate_Dn(i, Ai0, Ai)
            if self.continuity == 'C2':
                Di = [self.A0[i]] + list(Ai0) + [An, self.zetaT[i]]
            elif self.continuity == 'C1':
                Di = list(Ai0) + [An, self.zetaT[i]]
            # print(Di)
            self.cst[i].D = Di
            self.cst[i].zetaT = self.zetaT[i]
            self.cst[i].zetaL = self.zetaL[i]
            self.cst[i].offset_x = offset_x
            self.cst[i].internal_variables(self.cst[i].length)
            if i == 0:
                error = 0
            else:
                current = np.array([self.cst[i].chord, self.cst[i].zetaT,
                                    self.cst[i].zetaL, self.cst[i].D[0]])
                error = np.linalg.norm(current-prev)
                prev = current

    def _update_dependent(self, i):
        j = i - 1
        Ai0, Ai = self._select_D(i)

        error = 999
        k = 0
        while error > 1e-8:
            prev = np.array([self.cst[i].chord, self.cst[i].zetaT,
                             self.cst[i].zetaL, self.cst[i].D[0]])
            # print('p', i, k, prev)
            if i == 0:
                offset_x = 0
                self.A1[i] = self.D[0]
                self.cst[i].chord = self.spar_x[self.spar_i]
                self.zetaT[i] = (self.spar_y[self.spar_i]-self.offset)/self.cst[i].chord
                self.zetaL[i] = 0
                if self.root_fixed:
                    self.A0[i] = - self.zetaT[i]
                else:
                    self.A0[i] = self.D[0]
            elif self.N1[i] == 1. and self.N2[i] == 1:
                offset_x = self.cst[j].offset_x + self.cst[j].chord
                self.cst[i].chord = self.spar_x[self.spar_i] - offset_x
                self.zetaT[i] = (self.spar_y[self.spar_i]-self.offset)/self.cst[i].chord
                if self.continuity == 'C2':
                    ddj = self.n*self.cst[j].D[-3] - (self.N1[j]+self.n)*self.cst[j].D[-2]

                    self.A1[i] = (self.cst[i].chord/self.cst[j].chord *
                                  ddj+self.A0[i]*(self.n+1))/self.n
                elif self.continuity == 'C1':
                    raise(NotImplementedError)

                self.zetaL[i] = self.cst[j].chord/self.cst[i].chord*self.zetaT[j]
                self.A0[i] = -self.cst[j].D[-2] + self.zetaT[j] - \
                    self.zetaT[i] + self.zetaL[i] - self.zetaL[j]
            else:
                raise(NotImplementedError)
            self.cst[i]._D[0] = self.A0[i]
            self.cst[i]._D[1] = self.A1[i]
            self.cst[i]._D[-1] = self.zetaT[i]
            self.cst[i].zetaT = self.zetaT[i]
            self.cst[i].zetaL = self.zetaL[i]
            self.cst[i].offset_x = offset_x
            # print('outside 1', self.cst[i].D)
            An = self._calculate_Dn(i, Ai0, Ai)
            if self.continuity == 'C2':
                Di = [self.A0[i], self.A1[i]] + list(Ai0) + [An, self.zetaT[i]]
            elif self.continuity == 'C1':
                raise(NotImplementedError)

            self.cst[i].D = Di
            # print('before', i, self.cst[i].chord)
            # print('outside 2', self.cst[i].D)
            self.cst[i].internal_variables(self.cst[i].length)
            # print('after', i, self.cst[i].chord)
            if i == 0:
                error = 0
            else:
                # current = np.array([chord0, self.zetaT[i],
                #                     self.zetaL[i], self.A0[i]])
                current = np.array([self.cst[i].chord, self.cst[i].zetaT,
                                    self.cst[i].zetaL, self.cst[i].D[0]])
                error = np.linalg.norm(current-prev)

                # print('c', i, k, current)
                # print('error', i, k, error)
                prev = current
            k = k + 1
        self.spar_i += 1

    def _select_D(self, i):
        if i == 0:
            if self.continuity == 'C2':
                self.n_end = 1
            else:
                self.n_end = 0
            if self.dependent[i] and not self.root_fixed:
                self.n_end += 1
        if i == 1 and self.rigid_LE:
            self.n_end = 0
        self.n_start = self.n_end
        if self.dependent[i]:
            if self.length_preserving:
                modifier = 2
            else:
                modifier = 1
        else:
            modifier = 0
        if i == self.p-1 and self.free_end:
            self.n_end = self.n_start + self.nn - 1
            Ai0 = self.D[self.n_start:self.n_end]
            Ai = self.D[self.n_start:self.n_end]

        else:
            self.n_end = self.n_start + self.nn - modifier
            Ai = self.D[self.n_start:self.n_end]
            Ai0 = Ai[:-1]
        # print(i, self.n_start, self.n_end, Ai)
        return (Ai0, Ai)

    def _calculate_Dn(self, i, Ai0, Ai=None):
        def f(An):
            self.cst[i].D[-2] = An[0]
            self.cst[i].chord = chord_target
            self.cst[i].internal_variables(self.cst[i].length, origin=0)
            # length_current = self.cst[i].arclength(self.cst[i].chord)
            # print('Dn', self.cst[i].D, self.cst[i].chord, chord_target)
            return self.cst[i].chord - chord_target

        def fprime(An):
            self.cst[i]._D[-2] = An[0]
            return [integrate.quad(integrand, 0, chord_target, limit=500)[0]]

        def integrand(x):
            psi = x/self.cst[i].chord
            d = self.cst[i].x3(np.array([x]), diff='x1', offset=False)
            n = self.n
            dA = (1+psi)*(-psi*(2+n)+(1+n))
            return d*dA/(1+(d)**2)**(1/2)

        if self.dependent[i]:
            if self.length_preserving:
                # print(i)
                chord_target = self.cst[i].chord
                # r = np.linspace(-5, 5, 100)
                # rf = []
                # for j in range(len(r)):
                #     rf.append(f([r[j]]))
                # plt.figure()
                # plt.plot(r, rf)
                # plt.show()

                An = optimize.fsolve(f, self.cst[i].D[-2], fprime=fprime)[0]
            else:
                An = Ai[-1]
            # length_current = self.cst[i].arclength(self.cst[i].chord)
            # print('inside 1', self.cst[i].D, chord_target, self.cst[i].chord)
        elif i == self.p-1 and self.free_end:
            Pi = self.cst_p[i].D[:-1]
            chordp = self.cst_p[i].chord
            den_p = (1+(-Pi[-1] + self.cst_p[i].zetaT - self.cst_p[i].zetaL)**2)**(1.5)
            den_c = (1+(-self.cst[i].D[-2] + self.cst[i].zetaT -
                        self.cst[i].zetaL)**2)**(1.5)
            rho_p = (1/chordp)*(self.n*Pi[-2] - (self.N1[i]+self.n)*Pi[-1])/den_p
            An = (-den_c*self.cst[i].chord*rho_p + Ai0[-1]*self.n)/(self.n+1)
        else:

            An = Ai[-1]
        return An

    def _calculate_spar(self, i, debug=False):
        # Bernstein Polynomial order
        # dependent child
        dc = self.cst[i]
        # dependent parent
        dp = self.cst_p[i]
        # independent child
        ic = self.g_independent.cst[i]
        # independent parent
        ip = self.g_independent.cst_p[i]

        # Defining local variables
        N1 = dc.N1
        N2 = dc.N2
        A_ic = ic.D[0:-1]
        A_ip = ip.D[:-1]
        A_dp = dp.D[:-1]
        c_ip = ip.chord
        c_dp = dp.chord
        zT_ip = ip.D[-1]*c_ip
        zT_dp = dp.D[-1]*c_dp
        zL_ip = ip.zetaL*c_ip
        zL_dp = dp.zetaL*c_dp

        xi_upper_children = []
        ic.internal_variables(ic.length)
        c_ic = ic.chord
        zT_ic = ic.D[-1]*c_ic
        zL_ic = ic.zetaL*c_ic
        if debug:
            psi_spar = np.array([0.199984/c_ic])
            print('x/psi', psi_spar*c_ic, psi_spar)
        else:
            psi_spar = np.array([1])
        if i == 0:
            offset_x = 0
        else:
            offset_x = self.cst[i-1].offset_x + self.cst[i-1].chord
        # psi_upper_children = calculate_psi_goal(psi_spar, A_ip, A_ic, zT_ip,
        #                                         c_ip, c_ic, N1=N1, N2=N2,
        #                                         deltaz_goal=zT_ic,
        #                                         deltaL_baseline=zL_ip,
        #                                         deltaL_goal=zL_ic)

        # Calculate xi for upper children. Do not care about lower so
        # just gave it random shape coefficients
        xi_upper_children = CST(psi_spar, 1., deltasz=[
            zT_ic/c_ic, -zT_ic/c_ic], Al=A_ic, Au=A_ic, N1=N1, N2=N2,
            deltasLE=[zL_ic/c_ic, zL_ic/c_ic])
        xi_upper_children = xi_upper_children['u']

        self.spar_psi_upper = psi_spar
        self.spar_xi_upper = np.array(xi_upper_children) + ic.offset/c_ic
        # print('A', A_dp, A_ip)
        xi_parent = CST(psi_spar, 1., deltasz=[
            zT_ip/c_ip, -zT_dp/c_dp], Al=A_dp, Au=A_ip, N1=N1, N2=N2,
            deltasLE=[zL_ip/c_ip, zL_ip/c_ip])
        self.delta_P = c_ip*(xi_parent['u']-xi_parent['l']) + ip.offset - \
            dp.offset

        # print('delta', self.delta_P)
        # Claculate orientation for children
        s_j = calculate_spar_direction(
            psi_spar, A_ip, A_ic, zT_ip, c_ic, N1=N1, N2=N2, deltaz_goal=zT_ic,
            deltaL_baseline=zL_ip, deltaL_goal=zL_ic)
        self.spar_directions = s_j
        print('spar_directions', self.spar_psi_upper*c_ic, self.spar_directions,
              self.g_independent.x3(np.array([self.spar_psi_upper*c_ic]), 'theta1'),
              self.g_independent.x1(np.array([self.spar_psi_upper*c_ic]), 'theta1'))
        spar_x = psi_spar*c_ic - self.delta_P*s_j[0] + offset_x
        spar_y = xi_upper_children*c_ic + ic.offset - self.delta_P*s_j[1]
        if i == 0:
            self.spar_x = spar_x
            self.spar_y = spar_y
        else:
            self.spar_x = np.append(self.spar_x, spar_x)
            self.spar_y = np.append(self.spar_y, spar_y)

        # print('coordinates (upper)', self.spar_psi_upper *
        #       c_ic + ic.offset_x, self.spar_xi_upper*c_ic)
        # print('coordinates (lower)', self.spar_x, self.spar_y)

    @classmethod
    def cylindrical(cls, D, chord=1, color='k', configuration=0):
        c = cls(D, color=color, chord=chord, n=1, configuration=configuration)
        c.x3 = c._x3_cylindrical
        c.name = 'cylindrical'
        return c

    def _x3_poly(self, x1, diff=None, D=None, offset=None):
        """ z2 (checked)"""
        if D is None:
            D = self.D
        if diff is None:
            return(D[5]*x1**5 + D[4]*x1**4 + D[3]*x1**3 + D[2]*x1**2 + D[1]*x1 + D[0])
        elif diff == 'x1':
            return(5*D[5]*x1**4 + 4*D[4]*x1**3 + 3*D[3]*x1**2 + 2*D[2]*x1 + D[1])
        elif diff == 'x11':
            return(20*D[5]*x1**3 + 12*D[4]*x1**2 + 6*D[3]*x1 + 2*D[2])
        elif diff == 'x111':
            return(60*D[5]*x1**2 + 24*D[4]*x1 + 6*D[3])
        elif diff == 'theta1':
            return self.x3(x1, 'x1')*self.x1(x1, 'theta1')
        elif diff == 'theta11':
            return self.x3(x1, 'x11')*self.x1(x1, 'theta1')**2 + \
                self.x3(x1, 'x1')*self.x1(x1, 'theta11')
        elif diff == 'theta3':
            return(self.a[2, :, 2])
            # return(np.ones(len(x1)))
        elif diff == 'theta31' or diff == 'theta13':
            return(self.x1(x1, 'theta11'))
            # return(np.zeros(len(x1)))
        else:
            return(np.zeros(len(x1)))

    def _x3_CST(self, x1, diff=None, offset=True):
        A = self.D[:-1]
        if offset:
            x1 = x1 - self.offset_x
        psi = x1 / self.chord
        if diff is None:
            return(self.offset + CST(x1, self.chord, deltasz=self.zetaT*self.chord, Au=A,
                                     N1=self.N1, N2=self.N2, deltasLE=self.zetaL*self.chord))
        elif diff == 'x1':
            # print('x3 inside', self.zetaT, self.zetaL)
            d = dxi_u(psi, A, self.zetaT, N1=self.N1, N2=self.N2, zetaL=self.zetaL)
            if x1[0] < self.tol and self.N1 == 1:
                d[0] = +A[0] + self.zetaT - self.zetaL
            if abs(x1[-1] - self.chord) < self.tol and self.N2 == 1:
                d[-1] = -A[-1] + self.zetaT - self.zetaL
            return d
        elif diff == 'x11':
            dd = (1/self.chord)*ddxi_u(psi, A, N1=self.N1, N2=self.N2)
            if x1[0] < self.tol and self.N1 == 1 and self.N2 == 1:
                n = len(A) - 1
                dd[0] = (1/self.chord)*(-2*(n+1)*A[0] + 2*A[1]*n)
            if abs(x1[-1] - self.chord) < self.tol and self.N2 == 1:
                n = len(A) - 1
                dd[-1] = (1/self.chord)*(2*n*A[-2] - 2*(self.N1+n)*A[-1])
            return dd
        elif diff == 'theta1':
            return self.x3(x1, 'x1', False)*self.x1(x1, 'theta1', False)
        elif diff == 'theta11':
            return self.x3(x1, 'x11', False)*self.x1(x1, 'theta1', False)**2 + \
                self.x3(x1, 'x1', False)*self.x1(x1, 'theta11', False)
        elif diff == 'theta3':
            return(np.zeros(len(x1)))

    def _x3_pCST(self, x1, diff=None):
        output = np.zeros(len(x1))
        if len(x1) == 1:
            c_min = 0
            for i in range(self.p):
                c_max = c_min + self.cst[i].chord
                if (x1[0] >= (c_min-self.tol)) & (x1[0] <= (c_max+self.tol)):
                    output[0] = self.cst[i].x3(x1, diff=diff)
                    break
                c_min = c_max
        else:
            for i in range(self.p):
                output[self.cst[i].indexes] = self.cst[i].x3(x1[self.cst[i].indexes], diff=diff)
        return output

    def _x3_cylindrical(self, x1, diff=None, R=None):
        """ z2 (checked)"""
        if R is None:
            R = self.D[0]
        if diff is None:
            return(R - np.sqrt(R**2 - x1**2))
        elif diff == 'x1':
            return(x1/np.sqrt(R**2 - x1**2))
        elif diff == 'x11':
            return(R**2/np.sqrt(R**2 - x1**2)**3)
        elif diff == 'theta3':
            return(np.zeros(len(x1)))

    def r(self, x1=None, diff=None):
        if x1 is None:
            x1 = self.x1_grid
        else:
            if type(x1) == float:
                x1 = np.array([x1])
        if diff == 'theta1':
            output = np.array([self.x1(x1, 'x1', offset=False),
                               self.x3(x1, 'x1', offset=True)]).T
            output = np.einsum('ij,i->ij', output, self.x1(x1, 'theta1', offset=True))
        else:
            output = np.array([self.x1(x1, diff, offset=False),
                               self.x3(x1, diff, offset=True)]).T
            self.position = output
        return (output)

    def basis(self, x1=None, diff=None):
        if x1 is None:
            x1 = self.x1_grid

        if diff is None:
            self.a = np.zeros([3, len(x1), 3])
            self.a[0, :, :] = np.array([self.x1(x1, 'x1')*self.x1(x1, 'theta1'),
                                        [0]*len(x1),
                                        self.x3(x1, 'x1')*self.x1(x1, 'theta1')]).T
            # self.a[0,:,:] = np.array([[1]*len(x1),
            #                          [0]*len(x1),
            #                           [0]*len(x1)]).T
            self.a[1, :, :] = np.array([[0, 1, 0], ]*len(x1))
            self.a[2, :, :] = np.cross(self.a[0, :, :], self.a[1, :, :])

        elif diff == 'theta':
            # Most components are null
            self.da = np.zeros([3, 3, len(x1), 3])
            #  a1 diff theta1
            for i in range(3):
                for j in range(3):

                    # cross terms are null for 2D case
                    x11_1 = np.einsum('ij,i->ij', self.r(x1, 'x11'), self.x1(x1,
                                                                             'theta%d' % (i+1))*self.x1(x1, 'theta%d' % (j+1)))
                    x11_2 = np.einsum('ij,i->ij', self.r(x1, 'x1'),
                                      self.x1(x1, 'theta%d%d' % (i+1, j+1)))
                    x22_1 = np.einsum('ij,i->ij', self.r(x1, 'x22'), self.x2(x1,
                                                                             'theta%d' % (i+1))*self.x1(x1, 'theta%d' % (j+1)))
                    x22_2 = np.einsum('ij,i->ij', self.r(x1, 'x2'),
                                      self.x2(x1, 'theta%d%d' % (i+1, j+1)))
                    self.da[i, j, :, :] = x11_1+x11_2 + x22_1+x22_2

    def christoffel(self, i, j, k, order=1):
        if order == 1:
            gik_j = self.dA[i, k, j]
            gjk_i = self.dA[j, k, i]
            gij_k = self.dA[i, j, k]
            return .5*(gik_j + gjk_i - gij_k)
        elif order == 2:
            raise NotImplementedError

    def metric_tensor(self, diff=None):

        if diff is None:
            self.A = np.zeros([3, 3, len(self.x1_grid)])
            for i in range(3):
                for j in range(3):
                    self.A[i, j] = np.einsum('ij,ij->i', self.a[i, :], self.a[j, :])

        elif diff == 'theta':
            self.dA = np.zeros([3, 3, 3, len(self.x1_grid)])
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        self.dA[i, j, k] = np.einsum('ij,ij->i', self.da[i, k, :],
                                                     self.a[j, :]) + \
                            np.einsum('ij,ij->i', self.a[i, :],
                                      self.da[j, k, :])

    def curvature_tensor(self):
        self.B = np.zeros([2, 2, len(self.x1_grid)])
        for alpha in range(2):
            for beta in range(2):
                self.B[alpha, beta] = self.christoffel(alpha, beta, 2)

    def arclength(self, chord=None, origin=None):
        def integrand(x1):
            dr = self.x3(np.array([x1]), 'x1', offset=False)
            return np.sqrt(1 + dr**2)

        if chord is None:
            chord = self.chord
        if origin is None:
            origin = self.origin
        if chord == 0:
            return [0]
        return integrate.quad(integrand, origin, chord, limit=500)[0]

    def arclength_index(self, index):
        x1 = self.x1_grid[index]
        dr = self.x3(np.array([x1]), 'x1')
        if np.isnan(dr):
            if x1 == 0:
                dr = self.x3(np.array([self.tol]), 'x1')
            else:
                dr = self.x3(np.array([x1-self.tol]), 'x1')
        self.darc[index] = np.sqrt(1 + dr[0]**2)
        return integrate.trapz(self.darc[:index+1], self.x1_grid[:index+1]-self.offset_x)

    def improper_arclength_index(self, index):
        x1 = self.x1_grid[index]
        self.darc[index] = self.bounded_dr(x1)
        bounded_integral = integrate.trapz(self.darc[:index+1], self.x1_grid[:index+1])
        return bounded_integral + self.unbounded_integral(x1, 0)

    def arclength_chord(self):
        self.darc = np.ones(len(self.x1_grid))
        dr = self.x3(np.array([1e-7]), 'x1', offset=False)
        self.darc[0] = np.sqrt(1 + dr[0]**2)
        for index in range(1, len(self.x1_grid)-1):
            x1 = self.x1_grid[index]
            dr = self.x3(np.array([x1]), 'x1', offset=False)
            self.darc[index] = np.sqrt(1 + dr[0]**2)
        self.darc[-1] = np.sqrt(1 + (-self.D[-2]+self.D[-1]-self.zetaL)**2)
        return integrate.trapz(self.darc, self.x1_grid)

    def improper_arclength_chord(self):
        self.darc = np.zeros(len(self.x1_grid))
        for index in range(1, len(self.x1_grid)):
            x1 = self.x1_grid[index]
            self.darc[index] = self.bounded_dr(x1)
        bounded_integral = integrate.trapz(self.darc, self.x1_grid)
        return bounded_integral + self.unbounded_integral(self.chord)

    def bounded_dr(self, x):
        A = self.N1**2*self.D[0]**2
        n = self.n - 2
        if x == 0:
            return self.D[0]**2*(-n-1) + 1.5*n*self.D[0]*self.D[1]
        elif x == 1:
            # - np.sqrt(1 + B/np.sqrt(x))
            return np.sqrt(1 + (-self.D[-2]+self.D[-1])**2) - np.sqrt(1 + A/x)
        else:
            dr = self.x3(np.array([x]), 'x1', offset=False)
            return np.sqrt(1 + dr[0]**2) - np.sqrt(1 + A/x)

    def unbounded_integral(self, end, start=0):
        def indefinite_integral(x):
            return np.sqrt(x*(A + x)) + A*np.log(np.sqrt(A+x) + np.sqrt(x))
        A = self.N1**2*self.D[0]**2
        return indefinite_integral(end) - indefinite_integral(start)

    def calculate_x1(self, length_target, bounds=None, output=False, origin=0, length_rigid=0):
        def f(c_c):
            length_current = self.arclength(c_c[0])
            # print(c_c, length_current)
            return abs(target - length_current)

        def f_index(x):
            # Penalize in case x goes negative
            if x < 0:
                return 100
            else:
                self.x1_grid[index] = x + self.offset_x
                # if self.name == 'CST':
                #     length_current = length_rigid + self.improper_arclength_index(index)
                # else:
                length_current = length_rigid + self.arclength_index(index)
                return abs(target - length_current)

        if self.name == 'pCST':
            self.s = length_target
            rigid_n = 0
            for i in range(self.p):
                if i == 0 and self.rigid_LE:
                    rigid_n = len(self.cst[i].indexes)
                    self.cst[i].x1_grid = np.linspace(0, self.cst[i].chord,
                                                      rigid_n)
                else:
                    indexes = [self.cst[i].indexes[j] -
                               rigid_n for j in range(len(self.cst[i].indexes))]
                    print('length_target', length_target[indexes])
                    self.cst[i].s = length_target[indexes]
                    self.cst[i].calculate_x1(
                        length_target[indexes] - self.cst[i].offset_s)
            x1_grid = []
            for i in range(self.p):
                x1_grid += list(self.cst[i].x1_grid[:])
            if output:
                return x1_grid
            else:
                self.x1_grid = x1_grid
        else:
            x1 = []
            if len(length_target) == 1:
                target = length_target[0]
                x1.append(optimize.fsolve(f, origin)[0])
            else:
                self.x1_grid = np.zeros(len(length_target))
                self.darc = np.zeros(len(length_target))
                for index in range(len(length_target)):
                    if index == 0 and origin != 0:
                        dr = self.x3(np.array([origin + self.offset_x]), 'x1')
                        self.darc[index] = np.sqrt(1 + dr[0]**2)
                        self.x1_grid[index] = origin + self.offset_x
                    else:
                        target = length_target[index]
                        self.x1_grid[index] = optimize.fsolve(f_index, target)[0] + self.offset_x
            if output:
                return np.array(x1)

    def calculate_s(self, N, target_length=None, density='gradient', origin=0, output=False):
        def integrand(s):
            rho = self.radius_curvature(np.array([s]), output_only=True)[0]
            return abs(rho)

        def f(dx):
            if dx[0] >= 0:
                self.x1_grid[i] = self.x1_grid[i-1] + dx[0]
                self.rho[i] = integrand(self.x1_grid[i])
                partial = integrate.quad(integrand,
                                         self.x1_grid[i-1], self.x1_grid[i-1]+dx[0], limit=500)[0]
                self.partial[i] = partial
                return abs(partial - total/(N-1))
            else:
                return 100

        def fprime(dx):
            x = self.x1_grid[i-1]+dx[0]
            output = integrand(x)
            return [output]

        if self.name is not 'pCST' and target_length is None:
            target_length = self.length
        if density == 'gradient':
            if self.name == 'pCST':
                Nj = 0
                self.s = np.array([])
                self.indexes = []
                for i in range(self.p):
                    if i != 0 or not self.rigid_LE:
                        self.s = np.append(self.s, self.cst[i].calculate_s(N[i]))
                    self.cst[i].indexes = list(range(Nj, Nj + N[i]))
                    self.indexes += list(range(Nj, Nj + N[i]))
                    Nj += N[i]
            else:
                s_epsilon = self.arclength(np.array([origin]))[0] + self.offset_s
                self.s = np.linspace(s_epsilon, target_length + self.offset_s, N)
                return self.s
        elif density == 'curvature':
            total = integrate.quad(integrand, origin, self.chord, limit=500)[0]

            self.x1_grid = np.zeros(N)
            self.darc = np.zeros(N)
            self.partial = np.zeros(N)
            self.rho = np.zeros(N)
            self.rho[0] = abs(self.radius_curvature(np.array([0]), output_only=True)[0])
            s_list = [0]

            for i in range(1, N):
                dx = optimize.fsolve(f, 0, fprime=fprime)[0]
                self.x1_grid[i] = self.x1_grid[i-1] + dx
                s_list.append(self.improper_arclength_index(i))
                # BREAK
            return np.array(s_list)

    def plot(self, basis=False, r=None, label=None, linestyle='-', color=None, scatter=False, zorder=0, marker='.'):
        if self.name == 'pCST':
            for i in range(self.p):
                self.cst[i].plot(basis, label=label[i], scatter=scatter)
        else:
            if r is None:
                r = self.r(self.x1_grid)
            print('x_p', r[:, 0])
            print('y_p', r[:, 1])
            if color is None:
                color = self.color

            if scatter:
                if label is None:
                    plt.scatter(r[:, 0], r[:, 1], c=color, zorder=2, marker=marker)
                else:
                    plt.scatter(r[:, 0], r[:, 1], c=color, label=label,
                                zorder=2, edgecolors='k', marker=marker)
            else:
                if label is None:
                    plt.plot(r[:, 0], r[:, 1], color, linestyle=linestyle, lw=3,
                             zorder=zorder)
                else:
                    plt.plot(r[:, 0], r[:, 1], color, linestyle=linestyle, lw=3,
                             label=label, zorder=zorder)
            if basis:
                plt.quiver(r[:, 0], r[:, 2],
                           self.a[0, :, 0], self.a[0, :, 2],
                           angles='xy', color=color, scale_units='xy')
                plt.quiver(r[:, 0], r[:, 2],
                           self.a[2, :, 0], self.a[2, :, 2],
                           angles='xy', color=color, scale_units='xy')
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')

    def radius_curvature(self, x, output_only=False, parametric=False):
        if parametric:
            if self.name == 'CST' and x[0] == 0:
                x[0] = 1e-7
            dx = self.x1(x, diff='theta1')
            dy = self.x3(x, diff='theta1')
            ddx = self.x1(x, diff='theta11')
            ddy = self.x3(x, diff='theta11')
            rho = (dx*ddy-dy*ddx)/(dx**2 + dy**2)**(1.5)
        else:
            rho = self.x3(x, diff='x11')/(1+(self.x3(x, diff='x1'))**2)**(3/2)
            # if self.name == 'CST':
            #     if x[0] == 0:
            #         if self.D[0] == 0:
            #             rho[0] = 0
            #         else:
            #             rho[0] = -2/(self.D[0]**2)/self.chord
        if output_only:
            return rho
        else:
            self.rho = rho

    def internal_variables(self, target_length, origin=0):
        # At first we utilize the non-dimensional trailing edge thickness
        # print('inside 2', self.D)
        self.zetaT = self.D[-1]
        origin = origin/self.chord
        self.chord = 1

        nondimensional_length = self.arclength(chord=1., origin=origin/self.chord)

        self.chord = target_length/nondimensional_length
        self.deltaz = self.zetaT*self.chord
        self.deltazLE = self.zetaL*self.chord

    def calculate_angles(self):
        self.cos = self.x1(self.x1_grid, 'theta1')
        self.sin = self.x3(self.x1_grid, 'theta1')

        if self.x1_grid[0] == 0 and self.N1 == 1:
            self.cos[0] = 1
            self.sin[0] = 0

    def _check_input(self, input):
        # shape coefficients, zetaL, and zetaT
        total = (self.n+3)*self.p
        # y, dy boundary conditions
        dependent = 2*(self.p-1)
        # zetaL0 is always the same
        dependent += 1
        if self.free_end:
            dependent += 1
        if self.rigid_LE:
            dependent += self.n + 2
        else:
            if self.root_fixed:
                dependent += 1
        if self.continuity == 'C2':
            dependent += (self.p-1)
        # for structurally consistent

        if self.length_preserving:
            dependent += 2*np.count_nonzero(self.dependent)
        else:
            dependent += 1*np.count_nonzero(self.dependent)
        # print('Trues', np.count_nonzero(self.dependent))
        independent = total - dependent
        return len(input) == independent, len(input), independent
