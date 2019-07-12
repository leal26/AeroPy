import matplotlib.pyplot as plt
import numpy as np
import pickle

f = open('loudness.p', 'rb')
loudness = pickle.load(f)
f.close()

# f = open('output.p', 'rb')
# all_output = pickle.load(f)
# f.close()

f = open('outputs_prestrain_bump.p', 'rb')  #
data = pickle.load(f, encoding='latin1')
f.close()

displacements = {}
temperatures = {}
steps = ['Step-2', 'Step-3']
mid_index = int(len(data['U'][steps[0]][0])/2)
print(mid_index)
U0 = np.linalg.norm(data['U'][steps[0]][0][mid_index])
for step in steps:
    displacements[step] = []
    temperatures[step] = []

    for i in range(len(data['U'][step])):
        displacements[step].append(np.linalg.norm(data['U'][step][i][mid_index]) - U0)
        temperatures[step].append(data['NT11'][step][i][mid_index])


plt.figure()
plt.plot(data['Time']['Step-2'],
         loudness['Step-2'], 'b', label='Cooling')
plt.plot(1.0 + np.array(data['Time']['Step-3']),
         loudness['Step-3'], 'r', label='Heating')
plt.legend()
plt.show()

plt.figure()
plt.plot(data['Time']['Step-2'],
         displacements['Step-2'], 'b', label='Cooling')
plt.plot(1.0 + np.array(data['Time']['Step-3']),
         displacements['Step-3'], 'r', label='Heating')
plt.legend()
plt.show()

plt.figure()
plt.plot(data['Time']['Step-2'],
         temperatures['Step-2'], 'b', label='Cooling')
plt.plot(1.0 + np.array(data['Time']['Step-3']),
         temperatures['Step-3'], 'r', label='Heating')
plt.legend()
plt.show()

plt.figure()
plt.plot(temperatures['Step-2'],
         displacements['Step-2'], 'b', label='Cooling')
plt.plot(temperatures['Step-3'],
         displacements['Step-3'], 'r', label='Heating')
plt.legend()
plt.show()


plt.figure()
plt.plot(temperatures['Step-2'],
         loudness['Step-2'], 'b', label='Cooling')
plt.plot(temperatures['Step-3'],
         loudness['Step-3'], 'r', label='Heating')
plt.legend()
plt.show()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(temperatures['Step-2'],
         displacements['Step-2'], 'b--', label='Cooling')
ax1.plot(temperatures['Step-3'],
         displacements['Step-3'], 'b', label='Heating')
ax2.plot(temperatures['Step-2'],
         loudness['Step-2'], 'k--', label='Cooling')
ax2.plot(temperatures['Step-3'],
         loudness['Step-3'], 'k', label='Heating')
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Displacements (m)', color='b')
ax2.set_ylabel('Loudness (PLdB)', color='k')
plt.show()
