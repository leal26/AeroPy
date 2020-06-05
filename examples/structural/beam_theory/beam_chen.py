import matplotlib.pyplot as plt
import numpy as np
import pickle

from aeropy.structural.beam import beam_chen
from aeropy.structural.stable_solution import properties, loads
from aeropy.geometry.parametric import CoordinateSystem



g = CoordinateSystem.polynomial(D=[0, 0, 0, 0], chord = 1, color = 'b')
p = properties()
l = loads()
s = np.linspace(0, 1, 10)

b = beam_chen(g, p, l, s)
b.s_to_x()
b.find_deflection()

# Results from Abaqus
abaqus_data = pickle.load(open('neutral_line.p', 'rb'))

# Euler-Bernoulle
EB_solution = l.concentrated_load[0][-1]/(6*p.young*p.inertia) * \
    np.array([0, 0, 3, -1])

curve_EB = CoordinateSystem.polynomial(D=EB_solution, chord = 1, color  ='0.5')

curve_EB.calculate_x1(b.s)
eulerBernoulle = curve_EB.r(b.s)

[x,y,z] = eulerBernoulle.T
plt.plot(x,z, '.5', lw = 3, label='Euler-Bernoulle', linestyle = '-', zorder=0)
plt.plot(b.x, b.y, 'b', label='Chen', linestyle = '--', lw = 3)
plt.scatter(abaqus_data['coord'][0:401:40,0], abaqus_data['coord'][0:401:40,1], c='g', label='FEA', edgecolors='k', zorder = 10)
plt.legend()
plt.show()
