from scipy.optimize import approx_fprime
import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem

from aeropy.CST_2D import dxi_u, dxi_l, ddxi_u, ddxi_l
from aeropy.geometry.airfoil import CST


g = CoordinateSystem.CST(D=[0.11397826, 0.10433884, 0.10241407, 0.10070566, 0.0836374,  0.11353368, 0], chord=1,
                         color='b', N1=.5, N2=1, deltaz=0)
epsilon = 0.01
d = np.zeros([2, 1])
d[0] = g.x3(np.array([epsilon]))[0]
d[1] = g.x3(np.array([epsilon]), diff='x1')[0]

A0 = np.zeros(len(g.D)-1)
D = np.zeros([2, 2])
for i in range(2):
    A = np.copy(A0)
    A[i] = 1
    D[0, i] = CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=1)
    D[1, i] = dxi_u(epsilon, A, 0, N1=0.5, N2=1)

for i in range(2, 6):
    A = np.copy(A0)
    A[i] = g.D[i]
    d[0] -= CST(epsilon, Au=A, deltasz=0, N1=0.5, N2=1, c=1)
    d[1] -= dxi_u(epsilon, A, 0, N1=0.5, N2=1)
print(np.linalg.solve(D, d))
