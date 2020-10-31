import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.linalg import cho_factor, cho_solve

from aeropy.CST_2D import K, ddxi_u
from aeropy.geometry.parametric import CoordinateSystem

n = 3

m = n + 3

B = np.zeros((m, n+1))

for i in range(n+1):
    Bi = np.zeros(m)
    a2 = -n**2-3*n-2
    a1 = (2+2*i)*(n+1)
    a0 = -i**2-i
    Ki = K(i, n)
    # a0
    for j in range(n+1-i):
        p = j
        sign = (-1)**p
        B[j+i, i] += Ki*a0*sign*K(p, n-i)
    # a1
    for j in range(1, n+2-i):
        q = j-1
        sign = (-1)**q
        B[j+i, i] += Ki*a1*sign*K(q, n-i)
    # a2
    for j in range(2, n+3-i):
        r = j - 2
        sign = (-1)**r
        B[j+i, i] += Ki*a2*sign*K(r, n-i)

# linear solution
input = np.array([-1, 1])
a = np.zeros(m)
# a[:len(input)] = input
a[1:len(input)+1] += - input
a[2:len(input)+2] += input

x = np.linalg.lstsq(B, a)[0]
print('lstsq', x)
At = np.transpose(B)
Ata_inv = np.linalg.inv(np.matmul(At, B))
A_inv = np.matmul(Ata_inv, At)
x = np.matmul(A_inv, a)
print('Hermite?', x)

# Verification
g = CoordinateSystem.CST(D=list(x)+[0], chord=1, color='b', N1=1, N2=1)
g.calculate_s(21)
g.calculate_x1(g.s)
x1 = g.x1_grid
dd = g.x3(x1, diff='x11')

p = CoordinateSystem.polynomial(D=input, chord=1, color='b')
dd_p = p.x3(x1)

plt.figure()
plt.plot(x1, dd, 'r')
plt.plot(x1, dd_p, '--b')
plt.ylim([-2, 2])

# DD analytical
B = np.zeros(m)
A = x
for i in range(n+1):
    Bi = np.zeros(m)
    a2 = -n**2-3*n-2
    a1 = (2+2*i)*(n+1)
    a0 = -i**2-i
    Ki = K(i, n)
    # a0
    for j in range(n+1-i):
        p = j
        sign = (-1)**p
        Bi[j+i] += a0*sign*K(p, n-i)
        # print('a0', j+i, a0*sign*K(p, n-i))
    # a1
    for j in range(1, n+2-i):
        q = j-1
        sign = (-1)**q
        Bi[j+i] += a1*sign*K(q, n-i)
        # print('a1',j+i, a1*sign*K(q, n-i))
    # a2
    for j in range(2, n+3-i):
        r = j - 2
        sign = (-1)**r
        Bi[j+i] += a2*sign*K(r, n-i)
        # print('a2',j+i, a2*sign*K(r, n-i))
    B += A[i]*Ki*Bi

dd_a = np.zeros(len(x1))
for i in range(m):
    dd_a += B[i]*(x1**i)/(x1**2-x1)

plt.figure()
plt.plot(x1, dd, 'r')
plt.plot(x1, dd_a, '--b')
plt.ylim([-2, 2])
plt.show()
