import matplotlib.pyplot as plt
import numpy as np
import math

from aeropy.geometry.airfoil import create_x, CST

plt.figure()
A = [2/math.pi, 2/math.pi, 2/math.pi]
chord = 2/math.pi
deltaz = -2/math.pi*chord

x = create_x(chord, distribution='polar')
for i in range(len(A)):
    A_i = [0,0,0]
    A_i[i] = A[i]
    y = CST(x, chord, deltaz, A_i)
    plt.plot(x,y,'--')
y = CST(x, chord, deltaz, A, N1=1, N2=1)
plt.plot(x,y,'k')
plt.axis('equal')
plt.grid()
plt.xlim([0,1])
plt.xlabel('$\psi$', fontsize=14)
plt.ylabel(r'$\xi$', fontsize=14)
plt.show()


plt.figure()
A = [1.]
N1 = [.5, 1., .5]
N2 = [1., 1., .5]
deltaz = 0.0
x = create_x(1, distribution='polar')
for i in range(len(N1)):
    N1_i = N1[i]
    N2_i = N2[i]
    y = CST(x, 1., deltaz, A, N1=N1_i, N2=N2_i)
    plt.plot(x,y,'-')
y = CST(x, 1., deltaz, A)
plt.axis('equal')
plt.xlim([0,1])
plt.grid()
plt.xlabel('$\psi$', fontsize=14)
plt.ylabel(r'$\xi$', fontsize=14)
plt.show()