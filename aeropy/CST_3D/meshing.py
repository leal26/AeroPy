import numpy as np

def uniform_mesh_generator(mesh):
    '''Default genrator for uniform meshes. If 
       len(mesh)==2, upper and lower have same
       mesh'''

    # Defining dimensions matrix
    if len(mesh) == 2:
        dimensions = {'upper': list(mesh),
                      'lower': list(mesh)}
    elif len(mesh) == 4:
        dimensions = {'upper': list(mesh[:2]),
                      'lower': list(mesh[2:])}
    # Calculating grid
    upper_x, upper_y = np.meshgrid(np.linspace(0,1,dimensions['upper'][0]), 
                                   np.linspace(0,1,dimensions['upper'][1]))
    lower_x, lower_y = np.meshgrid(np.linspace(0,1,dimensions['lower'][0]), 
                                   np.linspace(0,1,dimensions['lower'][1]))

    mesh = {'upper': np.concatenate([upper_x.reshape(1,np.size(upper_x)),
                                     upper_y.reshape(1,np.size(upper_y))]).T,
            'lower': np.concatenate([lower_x.reshape(1,np.size(lower_x)),
                                     lower_y.reshape(1,np.size(lower_y))]).T}

    return mesh, dimensions