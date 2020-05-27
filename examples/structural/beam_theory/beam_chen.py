from scipy.integrate import quad
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np

P = 1
E = 1
L = 1
I = 1

def find_chord(P):
    def _to_minimize(l):
        def _to_integrate(y):
            den = np.sqrt(1-P**2/E**2/I**2*(l*y-y**2/2)**2)
            if np.isnan(den):
                return 100
            else:
                return 1/den
        l = l[0]
        current_L = quad(_to_integrate, 0, l)[0]
        return abs(L-current_L)
        
    return minimize(_to_minimize, L, method = 'Nelder-Mead',).x[0]

def find_deflection(y, l, P):
    def _to_integrate(y):
        num = P/E/I*(l*y-y**2/2)
        den = np.sqrt(1-P**2/E**2/I**2*(l*y-y**2/2)**2)
        if np.isnan(den):
            return 100
        else:
            return num/den
    
    x = []
    for y_i in y:
        x_i = current_L = quad(_to_integrate, 0, y_i)[0]
        x.append(x_i)
    return x
    
def find_chord_normalized(K):
    def _to_minimize(l):
        def _to_integrate(y):
            den = np.sqrt(1-K**2*(l*y-y**2/2)**2)
            if np.isnan(den):
                return 100
            else:
                return 1/den
        l = l[0]
        current_L = quad(_to_integrate, 0, l)[0]
        return abs(L-current_L)
        
    return minimize(_to_minimize, L, method = 'Nelder-Mead',).x[0]

def find_deflection_normalized(y, l, K):
    def _to_integrate(y):
        num = K*(l*y-y**2/2)
        den = np.sqrt(1-K**2*(l*y-y**2/2)**2)
        # print(num,den,K**2*(l*y-y**2/2)**2)
        if np.isnan(den):
            return 100
        else:
            return num/den
    
    x = []
    for y_i in y:
        print(y_i)
        x_i = current_L = quad(_to_integrate, 0, y_i)[0]
        x.append(x_i)
    return x
chord = find_chord(P)
y = np.linspace(0,chord,10)
x = find_deflection(y, chord, P)

plt.figure()
plt.plot(x,y)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('x')
plt.ylabel('y')

# For validation

K = np.linspace(0,20,50)
x = []
y = []
for K_i in K:
    chord = find_chord(K_i)
    deflection = find_deflection([chord], chord, K_i)[0]
    x.append(deflection)
    y.append(chord)
    print(K_i, chord, deflection)
plt.figure()
plt.plot(K, x, label='$\delta/L$')
plt.plot(K, y, label='$l/L$')
plt.legend()
plt.ylim([0, 1])
plt.xlim([0, 20])
plt.show()
