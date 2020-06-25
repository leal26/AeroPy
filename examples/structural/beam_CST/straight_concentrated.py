import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem

def format_input(input):
    return input[2:4]

g = CoordinateSystem.CST(D=[0, 0, 0, 0], chord = 1, color = 'b')
p = properties()
l = loads(concentrated_load = [[0, -1]], load_s = [1])
s = np.linspace(0, 1, 10)

b = beam_chen(g, p, l, s)
b.parameterized_solver(format_input)

plt.figure()
plt.plot(b.x, b.M)
plt.show()
# Results from Abaqus
abaqus_data = pickle.load(open('neutral_line.p', 'rb'))

# Euler-Bernoulle
EB_solution = l.concentrated_load[0][-1]/(6*p.young*p.inertia) * \
    np.array([0, 0, 3, -1])
print('EB', EB_solution)
curve_EB = CoordinateSystem.polynomial(D=EB_solution, chord = 1, color  ='0.5')
curve_EB.calculate_x1(b.s)
eulerBernoulle = curve_EB.r(b.s)

g0 = CoordinateSystem.polynomial(D=[0,0,0,0], chord = 1, color  ='k')
g0.calculate_x1(b.s)
g0.plot(label='Parent')
[x,y] = eulerBernoulle.T
plt.scatter(x,y, c='.5', label='Euler-Bernoulle: %.3f N' % -l.concentrated_load[0][1], linestyle = '-', zorder=1, edgecolors='k')
plt.scatter(abaqus_data['coord'][0:401:40,0], abaqus_data['coord'][0:401:40,1], c='.5', label='FEA: %.3f N' % -l.concentrated_load[0][1], edgecolors='k', zorder = 2, marker="^")
plt.plot(b.x, b.y, '.5', label='Child: %.3f N' % -l.concentrated_load[0][1], lw = 3, zorder=0)

plt.legend()
plt.show()
