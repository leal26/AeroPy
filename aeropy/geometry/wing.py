import math 
import numpy as np

def rotate(ux, uy, uz, theta):
    # Check if vector are normalized
    norm = math.sqrt(ux**2+uy**2+uz*2)
    if norm != 1:
        ux = ux/norm
        uy = uy/norm
        uz = uz/norm
    R11 = math.cos(theta) + ux**2*(1-math.cos(theta))
    R12 = ux*uy*(1-math.cos(theta)) - uz*math.sin(theta)
    R13 = ux*uz*(1-math.cos(theta)) + uy*math.sin(theta)
    R21 = ux*uy*(1-math.cos(theta)) + uz*math.sin(theta)
    R22 = math.cos(theta) + uy**2*(1-math.cos(theta))
    R23 = uy*uz*(1-math.cos(theta)) - ux*math.sin(theta)
    R31 = ux*uz*(1-math.cos(theta)) - uy*math.sin(theta)
    R32 = uy*uz*(1-math.cos(theta)) + ux*math.sin(theta)
    R33 = math.cos(theta) + uz**2*(1-math.cos(theta))
    R = np.array([[R11, R12, R13],
                  [R21, R22, R23],
                  [R31, R32, R33]])
    return R