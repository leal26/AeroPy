import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.linalg import cho_factor, cho_solve

from aeropy.CST_2D import K, ddxi_u
from aeropy.geometry.parametric import CoordinateSystem
from aeropy.structural.stable_solution import properties

n = 4

m = n + 3

D0 = np.array([0.20000000000000004, 0.2500000000000002, 0.2583333333333329,
               0.22500000000000048, 0.1500000000000001, -0.20000000000000004])
g = CoordinateSystem.CST(D=list(D0), chord=1, color='b', N1=1, N2=1)
g_p = CoordinateSystem.CST(D=list(D0), chord=1, color='b', N1=1, N2=1)
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

At = np.transpose(B)
Ata_inv = np.linalg.inv(np.matmul(At, B))
A_inv = np.matmul(Ata_inv, At)

# linear solution
input = np.array([-1, 1])
p = properties()
input = input/(p.young*p.inertia)

a0 = np.zeros(m)
# a[:len(input)] = input
a0[1:len(input)+1] += - input
a0[2:len(input)+2] += input

error = 9999
a = a0
x0 = [0]*(n+1)
while error > 1e-3:
    x = np.matmul(A_inv, a)
    g.D = list(D0 + np.array(list(x) + [-x[0]]))
    g.internal_variables(g.length)
    a = np.zeros(m)
    for k in range(m):
        a[k] = a0[k]*(g.chord**i)/(g.chord**(n+2)-g_p.chord**(n+2))
    residual = np.linalg.norm(np.matmul(B, x) - a)
    delta_A = np.linalg.norm(x-x0)
    # error = min([delta_A, residual])
    print('x', x)
    print('multi', np.matmul(B, x))
    print('a', a)
    print('chord', g.chord, g_p.chord)
    error = residual/np.linalg.norm(a)

    print(error, x)
    BREAK


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
