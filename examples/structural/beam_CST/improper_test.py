from scipy import integrate, optimize
import numpy as np


def integral(N):
    x = np.linspace(0.5, 1, N)
    dy = bounded(x)
    # print('dy', dy)
    bounded_integral = integrate.trapz(dy, x)
    return bounded_integral + unbounded(1, 0.5)


def bounded(x):
    output = 1/np.sqrt(x-x**2) - 1/np.sqrt(1 - x)
    if x[-1] == 1:
        output[-1] = 0
    return output


def unbounded(end, start=0):
    def indefinite_integral(x):
        return -2*np.sqrt(1-x)
    return indefinite_integral(end) - indefinite_integral(start)


for N in [4, 8, 10, 100, 1000]:
    print(N, integral(N))
