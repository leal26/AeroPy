import math
import numpy as np
import matplotlib.pyplot as plt

def hicks_henne(x, z, alfa):
    z = np.array(z)
    for j in range(len(x)):
        x_j = x[j]
        f_1 = math.sqrt(x_j)*(1-x_j)/math.exp(15*x_j)
        f_2 = math.sqrt(x_j)*(1-x_j)/math.exp(15*x_j)
        f_3 = math.sqrt(x_j)*(1-x_j)/math.exp(15*x_j)
        f_4 = math.sqrt(x_j)*(1-x_j)/math.exp(15*x_j)
        f_5 = math.sqrt(x_j)*(1-x_j)/math.exp(10*x_j)
        for i in range(len(alfa)):
            b_i = math.sin(math.pi*x[j]**(math.log(0.5)/math.log(x_list[i])))**t[i]
            z[j] += alfa[i]*b_i
    return z

def B(x, k, i, t):
    if k == 0:
       return 1.0 if t[i] <= x < t[i+1] else 0.0
    if t[i+k] == t[i]:
       c1 = 0.0
    else:
       c1 = (x - t[i])/(t[i+k] - t[i]) * B(x, k-1, i, t)
    if t[i+k+1] == t[i+1]:
       c2 = 0.0
    else:
       c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * B(x, k-1, i+1, t)
    return c1 + c2

def bspline(x, t, c, k):
    n = len(t) - k - 1
    assert (n >= k+1) and (len(c) >= n)
    return sum(c[i] * B(x, k, i, t) for i in range(n))    

# Bezier Curve
k = 2
t = [0, 1, 2, 3, 4, 5, 6]
c = [-1, 2, 0, -1]
x = np.linspace(1.5, 4.5, 50)
y_bezier = []
for x_i in x:
    y_bezier.append(bspline(x_i, t, c ,k))
    
# Hicks-Henne
y_hh = hicks_henne(x, y_bezier, [.1], [.0000001],[4.])    
plt.plot(x,y_bezier,label='Bezier')
plt.plot(x,y_hh,label='Hickes-Henne')
plt.xlabel('x')
plt.ylabel('z')
plt.legend()
plt.show()