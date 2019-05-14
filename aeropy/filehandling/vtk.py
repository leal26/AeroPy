import numpy as np


def generate_surface(data, filename='panair'):
    '''
    Function to generate vtk files from a panair input mesh
    INPUT :
    - data is a list of networks which are 3D arrays with the dimensions being
    columns, rows, and coordinates)
    - 'filename' is a string to use in filenames.
    For example 'panair' will result in files called 'panair_network_1', etc.

    OUTPUT :
    The function will produce one or several files, one for each network,
    in the folder it's run from.
    '''
    from evtk.hl import gridToVTK
    # TODO: Currently not working when meshes seeds are not the same

    # X, Y, Z = data
    X = np.ascontiguousarray(data[:, :, 0])
    Y = np.ascontiguousarray(data[:, :, 1])
    Z = np.ascontiguousarray(data[:, :, 2])
    gridToVTK(filename+'_network', X[:, :, None], Y[:, :, None], Z[:, :, None],
              cellData={'test': Z[:, :, None]})
    # def _write_network(points_array, multiple_networks=False):
    #     n_columns = int(points_array.shape[0])
    #     n_rows = int(points_array.shape[1])

    #     X = np.zeros((n_rows, n_columns, 1))
    #     Y = np.zeros((n_rows, n_columns, 1))
    #     Z = np.zeros((n_rows, n_columns, 1))

    #     for i in range(n_columns):
    #         for j in range(n_rows):
    #             X[j, i, 0] = points_array[i, j, 0]
    #             Y[j, i, 0] = points_array[i, j, 1]
    #             Z[j, i, 0] = points_array[i, j, 2]

    #     if multiple_networks:
    #         gridToVTK(filename+'_network'+str(n+1), X, Y, Z)
    #     else:
    #         gridToVTK(filename+'_network', X, Y, Z)

    # if type(data) == dict:
    #     networks = len(data.keys())
    # else:
    #     networks = len(data)

    # if type(data) != dict:
    #     try:
    #         # check to see if list of networks or just a single one
    #         check = data[0][0][0][0]
    #         for n in range(networks):
    #             points_array = data[n]
    #         _write_network(points_array, multiple_networks=True)
    #     except:
    #         _write_network(data, multiple_networks=False)
    # else:
    #     for n in range(networks):
    #         points_array = data[list(data.keys())[n]]
    #         _write_network(points_array, multiple_networks=True)


def generate_points(data, filename):
    # TODO: Still cannot visualize points well
    from evtk.hl import pointsToVTK
    x, y, z = data.T
    # Sometimes np.arrays could have manipulated to no longer
    # be c-continuous os we have to impose it
    x = np.ascontiguousarray(x)
    y = np.ascontiguousarray(y)
    z = np.ascontiguousarray(z)
    pointsToVTK(filename, x, y, z)
