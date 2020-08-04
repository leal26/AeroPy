import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.optimize import fsolve

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem
from aeropy.CST_2D import calculate_c_baseline, dxi_u
from aeropy.geometry.airfoil import create_x, rotate, CST

# def con(input):
#     g.D = input
#     g.internal_variables(target_length)
#
#     dd = (2*n*g.D[-3] - 2*(g.N1+n)*g.D[-2])
#     d = -g.D[-2] + g.D[-1]
#     rho = (1./g.chord)*dd/(1+d**2)**(3/2)
#     print('rho', rho, rho_p, g.D[-1])
#     return(abs(rho-rho_p))


def con(input):
    g.D = format_input(input, g_p=g_p)
    g.internal_variables(target_length)

    lhs = -2/(g.D[0]**2*g.chord) + 2/(g_p.D[0]**2*g_p.chord)
    rhs = -1*g.chord/p.young/p.inertia
    return(lhs-rhs)


# def format_input(input, g=None, g_p=None):
#     g.zetaT = input[-1]
#
#     g.deltaz = g.zetaT*g.chord
#     error = 999
#     while error > 1e-3:
#         chord0 = np.copy(g.chord)
#         d = np.zeros([2, 1])
#         d[0] = g_p.x3(np.array([epsilon]))[0] - epsilon*g.zetaT
#         d[1] = g_p.x3(np.array([epsilon]), diff='x1')[0] - g.zetaT
#
#         A0 = np.zeros(len(g.D)-1)
#         D = np.zeros([2, 2])
#         for i in range(2):
#             A = np.copy(A0)
#             A[i] = 1
#             # print(i, A)
#             D[0, i] = CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=g.chord)
#             D[1, i] = dxi_u(epsilon/g.chord, A, 0, N1=0.5, N2=1)
#
#         for i in range(2, 6):
#
#             A = np.copy(A0)
#             A[i] = input[i-2]
#             # print(i, A)
#             d[0] -= CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=g.chord)
#             d[1] -= dxi_u(epsilon/g.chord, A, 0, N1=0.5, N2=1)
#         # A01 = np.linalg.solve(D, d)
#
#         det = D[0, 0]*D[1, 1] - D[0, 1]*D[1, 0]
#         A0 = (1/det)*(D[1, 1]*d[0][0] - D[0, 1]*d[1][0])
#         A1 = (1/det)*(-D[1, 0]*d[0][0] + D[0, 0]*d[1][0])
#         g.D = [A0, A1] + list(input)
#         g.internal_variables(b.length, origin=epsilon)
#         error = abs(g.chord-chord0)
#         print('iteration', error, g.chord, chord0)
#     return [A0, A1] + list(input)

def format_input(input, g=None, g_p=None):
    error = 999
    while error > 1e-3:
        chord0 = np.copy(g.chord)

        A0 = np.zeros(len(g.D)-1)
        # D
        D = np.zeros([2, 2])
        A = np.copy(A0)
        A[0] = 1
        D[0, 0] = CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=g.chord)
        D[1, 0] = dxi_u(epsilon/g.chord, A, 0, N1=0.5, N2=1)
        D[0, 1] = epsilon
        D[1, 1] = 1

        # d
        d = np.zeros([2, 1])
        d[0] = g_p.x3(np.array([epsilon]))[0]
        d[1] = g_p.x3(np.array([epsilon]), diff='x1')[0]
        for i in range(1, 6):
            A = np.copy(A0)
            A[i] = input[i-1]
            d[0] -= CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=g.chord)
            d[1] -= dxi_u(epsilon/g.chord, A, 0, N1=0.5, N2=1)

        det = D[0, 0]*D[1, 1] - D[0, 1]*D[1, 0]
        A0 = (1/det)*(D[1, 1]*d[0][0] - D[0, 1]*d[1][0])
        # print('D', D)
        # print('d', d)
        # BREAK
        b.g.zetaT = (1/det)*(-D[1, 0]*d[0][0] + D[0, 0]*d[1][0])
        g.D = [A0] + list(input) + [b.g.zetaT]
        g.internal_variables(b.length, origin=epsilon)
        error = abs(g.chord-chord0)
        print('iteration', error, g.chord, chord0)
    return [A0] + list(input) + [b.g.zetaT]

# abaqus_x = [0.0, 0.0099999998, 0.0020458214, 0.047780406, 0.087024547, 0.12658319, 0.16626321, 0.20599827, 0.24575868, 0.2855289, 0.32529992, 0.3650662, 0.40482411, 0.44457138,
#             0.48430675, 0.52402967, 0.56374025, 0.60343891, 0.64312649, 0.68280399, 0.72247189, 0.76213056, 0.80177915, 0.84141576, 0.88103628, 0.92063403, 0.96019793, 0.99971026]
# abaqus_x = [0.0, 0.0099999998, 0.00206689, 0.01915412, 0.028679075, 0.038376939, 0.04817313, 0.058031909, 0.067933477, 0.077865809, 0.087821044, 0.097793795, 0.10778025, 0.11777757, 0.12778363, 0.13779677, 0.14781572, 0.15783946, 0.16786712, 0.17789805, 0.18793166, 0.19796748, 0.20800516, 0.21804431, 0.22808467, 0.23812598, 0.24816802, 0.25821066, 0.26825362, 0.27829689, 0.28834024, 0.29838359, 0.30842686, 0.31846994, 0.32851276, 0.33855525, 0.34859732, 0.35863897, 0.36868006, 0.37872061, 0.3887606, 0.39879993, 0.40883863, 0.41887659, 0.42891383, 0.43895039, 0.44898614, 0.45902112, 0.46905535, 0.47908878,
#             0.48912144, 0.49915326, 0.5091843, 0.51921457, 0.52924401, 0.53927261, 0.54930055, 0.55932766, 0.569354, 0.57937956, 0.58940434, 0.59942853, 0.60945195, 0.61947465, 0.62949669, 0.63951808, 0.64953876, 0.65955889, 0.66957825, 0.67959714, 0.68961537, 0.69963294, 0.70964998, 0.71966648, 0.72968233, 0.73969758, 0.74971223, 0.75972629, 0.76973975, 0.77975249, 0.78976458, 0.79977596, 0.80978662, 0.8197965, 0.82980543, 0.83981353, 0.84982067, 0.85982662, 0.86983144, 0.87983501, 0.88983715, 0.89983767, 0.90983647, 0.91983342, 0.92982823, 0.93982065, 0.9498105, 0.95979744, 0.96978122, 0.97976136, 0.98973757, 0.99970925]
# abaqus_y = [0, 0.011236889, 0.00514034, 0.023266012, 0.029659046, 0.033741537, 0.036422059, 0.038113344, 0.039044935, 0.039360147, 0.039156213, 0.038503978, 0.037458673, 0.036066201,
#             0.034366865, 0.032397442, 0.03019212, 0.027782552, 0.025197186, 0.022459945, 0.019588264, 0.016590495, 0.013463141, 0.010186969, 0.0067228526, 0.0030072413, -0.0010527908, -0.0055848872]
# abaqus_y = [0, 0.011236889, 0.005166586, 0.015349808, 0.018527772, 0.021135807, 0.023348065, 0.025263296, 0.026943931, 0.028432477, 0.029759483, 0.030947819, 0.032015145, 0.032975465, 0.033840097, 0.034618367, 0.035318051, 0.035945732, 0.036507033, 0.037006792, 0.037449226, 0.037838027, 0.038176447, 0.038467374, 0.038713377, 0.038916763, 0.039079618, 0.039203819, 0.03929108, 0.039342966, 0.039360911, 0.039346233, 0.039300159, 0.039223831, 0.039118305, 0.038984574, 0.038823579, 0.038636211, 0.0384233, 0.038185652, 0.037924036, 0.037639186, 0.037331808, 0.037002593, 0.0366522, 0.036281284, 0.035890464, 0.035480365, 0.035051588, 0.034604717, 0.034140341, 0.033659026,
#             0.033161335, 0.032647807, 0.032119002, 0.031575438, 0.031017635, 0.030446114, 0.029861365, 0.029263875, 0.028654126, 0.028032573, 0.027399652, 0.026755802, 0.026101431, 0.02543691, 0.024762599, 0.024078859, 0.023386011, 0.022684297, 0.021973955, 0.021255219, 0.020528316, 0.019793352, 0.019050363, 0.018299377, 0.017540419, 0.016773501, 0.01599847, 0.015215059, 0.014423011, 0.013622062, 0.012811875, 0.011991893, 0.011161525, 0.010320176, 0.009467124, 0.0086014681, 0.0077222902, 0.006828553, 0.0059190569, 0.004992547, 0.004047587, 0.0030826547, 0.0020960681, 0.0010860161, 5.0512786e-05, -0.0010124723, -0.0021053595, -0.0032304232, -0.0043901685, -0.0055878009]
# rotated_abaqus = rotate({'x': abaqus_x, 'y': abaqus_y}, normalize=False)
# abaqus_x = rotated_abaqus['x']
# abaqus_y = rotated_abaqus['y']


# 1%
# abaqus_x = [0, 0.0046626413, 0.018402023, 0.048702486, 0.087927498, 0.12745513,
#             0.16709827, 0.20679522, 0.24651785, 0.28625035, 0.3259826,
#             0.36570814, 0.40542319, 0.44512612, 0.48481697, 0.52449685,
#             0.5641672, 0.60382956, 0.64348471, 0.68313259, 0.72277194,
#             0.76240051, 0.80201453, 0.84160942, 0.88117975, 0.92072034,
#             0.96022779, 0.99970293]
# abaqus_y = [0, 0.0080771167, 0.015637144, 0.024132574, 0.030406451,
#             0.034420185, 0.037082944, 0.038775206, 0.039692041, 0.039953224,
#             0.039647505, 0.038851682, 0.037638135, 0.036076538, 0.034232546,
#             0.032165118, 0.029923303, 0.027543247, 0.025045833, 0.022435322,
#             0.019699335, 0.016810443, 0.013729585, 0.010411576, 0.0068127746,
#             0.0029010926, -0.0013316774, -0.0058551501]
# abaqus_x = [0.0, 0.0099999998, 0.00206689, 0.01915412, 0.028679075, 0.038376939, 0.04817313, 0.058031909, 0.067933477, 0.077865809, 0.087821044, 0.097793795, 0.10778025, 0.11777757, 0.12778363, 0.13779677, 0.14781572, 0.15783946, 0.16786712, 0.17789805, 0.18793166, 0.19796748, 0.20800516, 0.21804431, 0.22808467, 0.23812598, 0.24816802, 0.25821066, 0.26825362, 0.27829689, 0.28834024, 0.29838359, 0.30842686, 0.31846994, 0.32851276, 0.33855525, 0.34859732, 0.35863897, 0.36868006, 0.37872061, 0.3887606, 0.39879993, 0.40883863, 0.41887659, 0.42891383, 0.43895039, 0.44898614, 0.45902112, 0.46905535, 0.47908878,
#             0.48912144, 0.49915326, 0.5091843, 0.51921457, 0.52924401, 0.53927261, 0.54930055, 0.55932766, 0.569354, 0.57937956, 0.58940434, 0.59942853, 0.60945195, 0.61947465, 0.62949669, 0.63951808, 0.64953876, 0.65955889, 0.66957825, 0.67959714, 0.68961537, 0.69963294, 0.70964998, 0.71966648, 0.72968233, 0.73969758, 0.74971223, 0.75972629, 0.76973975, 0.77975249, 0.78976458, 0.79977596, 0.80978662, 0.8197965, 0.82980543, 0.83981353, 0.84982067, 0.85982662, 0.86983144, 0.87983501, 0.88983715, 0.89983767, 0.90983647, 0.91983342, 0.92982823, 0.93982065, 0.9498105, 0.95979744, 0.96978122, 0.97976136, 0.98973757, 0.99970925]
# abaqus_y = [-4.0000065e-40, 0.011236889, 0.005166586, 0.015349808, 0.018527772, 0.021135807, 0.023348065, 0.025263296, 0.026943931, 0.028432477, 0.029759483, 0.030947819, 0.032015145, 0.032975465, 0.033840097, 0.034618367, 0.035318051, 0.035945732, 0.036507033, 0.037006792, 0.037449226, 0.037838027, 0.038176447, 0.038467374, 0.038713377, 0.038916763, 0.039079618, 0.039203819, 0.03929108, 0.039342966, 0.039360911, 0.039346233, 0.039300159, 0.039223831, 0.039118305, 0.038984574, 0.038823579, 0.038636211, 0.0384233, 0.038185652, 0.037924036, 0.037639186, 0.037331808, 0.037002593, 0.0366522, 0.036281284, 0.035890464, 0.035480365, 0.035051588, 0.034604717, 0.034140341,
#             0.033659026, 0.033161335, 0.032647807, 0.032119002, 0.031575438, 0.031017635, 0.030446114, 0.029861365, 0.029263875, 0.028654126, 0.028032573, 0.027399652, 0.026755802, 0.026101431, 0.02543691, 0.024762599, 0.024078859, 0.023386011, 0.022684297, 0.021973955, 0.021255219, 0.020528316, 0.019793352, 0.019050363, 0.018299377, 0.017540419, 0.016773501, 0.01599847, 0.015215059, 0.014423011, 0.013622062, 0.012811875, 0.011991893, 0.011161525, 0.010320176, 0.009467124, 0.0086014681, 0.0077222902, 0.006828553, 0.0059190569, 0.004992547, 0.004047587, 0.0030826547, 0.0020960681, 0.0010860161, 5.0512786e-05, -0.0010124723, -0.0021053595, -0.0032304232, -0.0043901685, -0.0055878009]
# 10%
abaqus_x = [0.0, 0.1, 0.0020865877, 0.010215517, 0.019589836, 0.029336579, 0.03925629, 0.049274135, 0.059354335, 0.069477111, 0.079630449, 0.089806534, 0.10996284, 0.11993612, 0.12991785, 0.13990647, 0.14990073, 0.15989965, 0.16990244, 0.17990845, 0.18991712, 0.19992799, 0.20994069, 0.21995488, 0.22997028, 0.23998664, 0.25000378, 0.26002148, 0.27003959, 0.280058, 0.29007655, 0.30009511, 0.31011364, 0.32013196, 0.3301501, 0.34016794, 0.35018536, 0.3602024, 0.37021899, 0.38023502, 0.39025047, 0.40026537, 0.41027957, 0.42029315, 0.43030608, 0.44031826, 0.45032972, 0.46034044, 0.47035041, 0.48035964, 0.4903681,
            0.50037575, 0.51038271, 0.52038878, 0.5303942, 0.54039878, 0.55040264, 0.56040579, 0.57040811, 0.58040977, 0.59041071, 0.60041094, 0.61041045, 0.62040931, 0.63040751, 0.64040506, 0.65040201, 0.6603983, 0.670394, 0.68038905, 0.69038355, 0.70037746, 0.71037078, 0.72036356, 0.7303558, 0.74034739, 0.75033838, 0.76032883, 0.77031863, 0.78030783, 0.79029632, 0.80028415, 0.8102712, 0.82025754, 0.83024305, 0.8402276, 0.8502112, 0.86019373, 0.87017512, 0.88015527, 0.89013398, 0.90011114, 0.91008663, 0.92006034, 0.93003184, 0.94000107, 0.9499678, 0.95993161, 0.96989226, 0.97984946, 0.98980272, 0.99975163]
abaqus_y = [-4.0000065e-40, 0.031264387, 0.0051910011, 0.011353852, 0.015514729, 0.018725043, 0.021357819, 0.023590373, 0.025523078, 0.027219331, 0.028722271, 0.030062895, 0.032319341, 0.033270512, 0.034128651, 0.03490264, 0.035599921, 0.036226809, 0.036788702, 0.037290271, 0.037735589, 0.038128216, 0.038471311, 0.038767669, 0.039019786, 0.039229915, 0.039400075, 0.039532099, 0.039627664, 0.039688289, 0.039715387, 0.039710242, 0.039674055, 0.039607946, 0.039512955, 0.039390061, 0.039240185, 0.039064191, 0.038862918, 0.038637146, 0.038387638, 0.038115114, 0.03782028, 0.037503809, 0.037166361, 0.036808569, 0.036431063, 0.036034454, 0.035619333, 0.035186291, 0.034735892,
            0.034268714, 0.033785302, 0.033286206, 0.03277196, 0.032243092, 0.031700108, 0.031143531, 0.030573845, 0.029991528, 0.029397061, 0.028790889, 0.02817345, 0.027545173, 0.026906464, 0.026257688, 0.025599198, 0.024931356, 0.024254471, 0.023568785, 0.022874529, 0.022171937, 0.021461232, 0.020742511, 0.020015802, 0.019281134, 0.018538537, 0.017788017, 0.017029405, 0.01626244, 0.015486865, 0.01470242, 0.013908756, 0.013105318, 0.012291517, 0.011466763, 0.010630327, 0.009781315, 0.0089188144, 0.008041786, 0.0071490384, 0.0062393229, 0.0053112092, 0.0043631843, 0.003393573, 0.0024005773, 0.0013822207, 0.00033648894, -0.00073902967, -0.0018465805, -0.0029886768, -0.0041684718]
epsilon = 0.1
g = CoordinateSystem.CST(D=[0.11397826, 0.10433884, 0.10241407, 0.10070566, 0.0836374, 0.11353368, 0.], chord=1,
                         color='b', N1=.5, N2=1)
g.name = 'proper integral'
g_p = CoordinateSystem.CST(D=[0.11397826, 0.10433884, 0.10241407, 0.10070566, 0.0836374, 0.11353368, 0], chord=1,
                           color='k', N1=.5, N2=1)
g_p.name = 'proper integral'
s_epsilon = g_p.arclength(np.array([epsilon]))[0]
s = np.linspace(s_epsilon, g_p.arclength(np.array([1]))[0], 51)
# s = g_p.calculate_s(11, density='curvature', origin=epsilon)
p = properties(dimensions=[0.001, 0.001])
l = loads(concentrated_load=[[0, -0.0001]], load_s=[s[-1]])
print('s', s)
b = beam_chen(g, p, l, s, origin=epsilon, ignore_ends=True)

n = g_p.n - 2
dd_p = (2*n*g_p.D[-3] - 2*(g_p.N1+n)*g_p.D[-2])
d_p = -g_p.D[-2] + g_p.D[-1]
rho_p = (1/g_p.chord)*dd_p/(1+d_p**2)**(3/2)
target_length = s[-1]

# cons = {'type': 'eq', 'fun': con}
# n = len(b.g_p.D) - 2
# Cn0 = g.D[-2] - 2*n/(1+2*n)*g.D[-3]
b.parameterized_solver(format_input=format_input, x0=g.D[1:-1])

print('x', b.x)
print('y', b.y)
print('r', b.r)
print('R', b.R)
print('arc', b.g.darc)
print('rho', b.g.rho)

b.g.calculate_x1(b.s, origin=b.origin, length_rigid=b.s[0])
b.x = b.g.x1_grid
b.y = b.g.x3(b.x)
b.g.radius_curvature(b.g.x1_grid)
b.g_p.radius_curvature(b.g_p.x1_grid)
# rotated_beam = rotate({'x': b.x, 'y': b.y}, normalize=False)
plt.figure()
plt.plot(b.x, b.M/b.p.young/b.p.inertia, c='b')
plt.plot(b.x, b.g.rho, c='r')
plt.plot(b.x, b.g_p.rho, c='g')
plt.figure()
plt.plot(b.x, b.M/b.p.young/b.p.inertia, c='b')
plt.plot(b.x, b.g.rho - b.g_p.rho, c='r')
plt.show()
print('residual: ', b._residual(b.g.D))
print(b.r)
print('x', b.x)
print('y', b.y)
print('length', target_length, b.length, b.g.arclength(b.g.chord, origin=epsilon)[0])
print('dys', b.g_p.x3(np.array([b.g_p.x1_grid[0]]), 'x1'), b.g.x3(np.array([b.g.x1_grid[0]]), 'x1'))
print('ddys', b.g.x3(np.array([b.g.x1_grid[-1]]), 'x11'),
      b.g.x3(np.array([b.g.x1_grid[-1]]), 'x11'))
b.g_p.plot(label='Parent')
plt.plot(b.x, b.y, '.5',
         label='Child: %.3f N' % -l.concentrated_load[0][-1], lw=3)
plt.scatter(abaqus_x, abaqus_y, c='.5', label='FEA: %.3f N' % -
            l.concentrated_load[0][-1], edgecolors='k', zorder=10, marker="^")
# plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.show()
