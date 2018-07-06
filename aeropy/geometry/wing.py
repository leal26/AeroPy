import math 
import numpy as np

def Rotation_euler_vector(euler_vector, theta):
    ux = euler_vector[0]
    uy = euler_vector[1]
    uz = euler_vector[2]
    # Check if vector are normalized
    norm = math.sqrt(ux**2+uy**2+uz*2)
    if norm != 1:
        ux = ux/norm
        uy = uy/norm
        uz = uz/norm

    R11 = np.cos(theta) + ux**2*(1-np.cos(theta))
    R12 = ux*uy*(1-np.cos(theta)) - uz*np.sin(theta)
    R13 = ux*uz*(1-np.cos(theta)) + uy*np.sin(theta)
    R21 = ux*uy*(1-np.cos(theta)) + uz*np.sin(theta)
    R22 = np.cos(theta) + uy**2*(1-np.cos(theta))
    R23 = uy*uz*(1-np.cos(theta)) - ux*np.sin(theta)
    R31 = ux*uz*(1-np.cos(theta)) - uy*np.sin(theta)
    R32 = uy*uz*(1-np.cos(theta)) + ux*np.sin(theta)
    R33 = np.cos(theta) + uz**2*(1-np.cos(theta))
    R = np.array([[R11, R12, R13],
                  [R21, R22, R23],
                  [R31, R32, R33]])
    return R

def Rotation_euler_angles(theta):
 
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                    ])
                     
    R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                    ])
                 
    R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])   
                     
    R = np.dot(R_z, np.dot( R_y, R_x ))
 
    return R