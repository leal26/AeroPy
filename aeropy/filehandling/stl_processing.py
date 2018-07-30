import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import stl
import pickle
import os

from aeropy.xfoil_module import output_reader
import aeropy.geometry.airfoil as am

def points_from_stl(filename) :
    '''Takes the path to the file as string input
    Returns an array of points'''
    mesh = stl.mesh.Mesh.from_file(filename)
    list_points = []
    counter = 0
    for p in mesh.v0 :
        list_points.append(list(p))
    for p in mesh.v1 :
        list_points.append(list(p))
    for p in mesh.v2 :
        list_points.append(list(p))
        
    return np.unique(list_points, axis=0)


def extract_stl(directory = "JAXA_files", filenames = ['wing_left',
                'wing_left_cap','wing_right', 'wing_right_cap',
                'fuselage','full_coarse'], return_output = False,
                normalize = True, store_output = True):
    """Open Stl files and store information in pickles or return
    variable. Normalize divides everything by the half_span."""
    if return_output:
        output = {}
    for name in filenames:
        data = points_from_stl(directory+"%s.stl" % name)
        if store_output and store_output:
            pickle.dump( data, open( directory+'Raw_'+name+'.p', "wb" ) )
        if normalize:
            norm = max(data[:,1])
            data /= norm
        if return_output:
            output[name] = data
    if store_output and normalize:
        pickle.dump( data, open( directory+'normalized_data.p', "wb" ) )
    if return_output:
        return(output)
 
def filter_geometry(data, parts ):
    '''Filter part points to only select y-wise regions with considerable
       number of points'''
    def _high_density_eta(data, min_repetitions=100):
        '''Find values of eta with high density to use for leading edge'''
        TOL = 1e-3
        a = data[:,1]
        b = a.copy()
        b.sort()
        d =  np.append(True, np.diff(b))
        eta_unique = b[d>TOL]
        
        # Find values that have multiple repetitions
        eta_dense = []
        for eta in eta_unique:
            c = np.abs(a - eta)
            r = c[c<=TOL]
            if len(r) >= min_repetitions:
                print(len(r))
                eta_dense.append(eta)

        return eta_dense
    
    data_sampled = {}
    eta_sampled = {}
    tol=1e-3
    for part in parts:
        # Getting airfoil data for certain spanwise locations with high density of points
        eta_sampled[part] = _high_density_eta(data['wing_right'], 2)
        data_sampled[part] = {}
        for eta in eta_sampled[part]:
            data_sampled[part][eta] = []
        for i in range(len(data[part])):
            diff = min(np.abs(np.array(eta_sampled[part])-data[part][i][1]))
            index_sampled = np.where(np.abs(np.array(eta_sampled[part])-data[part][i][1])==diff)[0][0]
            if diff < tol:
                data_sampled[part][eta_sampled[part][index_sampled]].append(data[part][i])
        # convert to numpy arrays
        for eta in eta_sampled[part]:
            data_sampled[part][eta] = np.array(data_sampled[part][eta])
    return data_sampled, eta_sampled
 


def clean_edges(eta_dense, LE, TE):
    '''Based on the dervative, remove points that make edges look bad'''
    new_LE = []
    new_TE = []
    # Flag to know if sweep is negative or positive
    if (LE[-1][0]-LE[0][0]) > 0:
        sweep_direction_LE = 1
    else:
        sweep_direction_LE = -1
    if (TE[-1][0]-TE[0][0]) > 0:
        sweep_direction_TE =1
    else:
        sweep_direction_TE = -1

    # Filter points that have derivative that should not be
    LE = np.array(LE)
    TE = np.array(TE)
    
    x,y,z = LE.T
    diff_LE = np.append(sweep_direction_LE, np.diff(x))
    eta_dense = np.array(eta_dense)[sweep_direction_LE*diff_LE>=0]
    LE = LE[sweep_direction_LE*diff_LE>=0]
    TE = TE[sweep_direction_LE*diff_LE>=0]
    
    x,y,z = TE.T
    diff_TE = np.append(sweep_direction_LE, np.diff(x))    
    eta_dense = np.array(eta_dense)[sweep_direction_TE*diff_TE>=0]
    LE = LE[sweep_direction_TE*diff_TE>=0]
    TE = TE[sweep_direction_TE*diff_TE>=0]

    return LE, TE, eta_dense

def edge_points(data_sampled, eta_sampled, parts):
    LE = {}
    TE = {}
    LE_list = []
    TE_list = []
    theta_list = []
    chord_list = []

    for part in parts:
        for eta in eta_sampled[part]:
            x,y,z = data_sampled[part][eta].T
            [LE_i,TE_i,theta,chord] = am.find_edges(list(x), list(z))
            LE_list.append([LE_i['x'], eta, LE_i['y'],])
            TE_list.append([TE_i['x'], eta, TE_i['y'],])

            theta_list.append(theta)
            chord_list.append(chord)
        [LE[part], TE[part], eta_sampled[part]] = clean_edges(eta_sampled[part], LE_list, TE_list) 

    pickle.dump([LE, TE], open( directory+"edges.p", "wb" ) )
    return LE, TE
    
if __name__ == '__main__': 
    directory = "./JAXA_files/"
    filenames = ['wing_right']
    extract_stl(directory)

    directory = "./examples/"
    filename = 'CST_twist'
    parts = ['wing_right']
    data = {}
    
    data_raw = pickle.load( open( directory + filename + ".p",  "rb" ))
    for part in parts:
        # Mixing lower and upper surface to test the LE and TE finder algorithms
        data[part] = np.append(data_raw['upper'], data_raw['lower'], axis=0)
        
    data_sampled, eta_sampled = filter_geometry(data,parts)
            
    LE, TE = edge_points(data_sampled, eta_sampled, parts)

    print(LE)
    # Printing just to check it out
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    for part in parts:
        for eta in eta_sampled[part]:
            print(data_sampled[part][eta])
            x,y,z = data_sampled[part][eta].T
            ax.scatter(x, z, label=r'$\eta$=%.3f' % eta)
        x,y,z=LE[part].T
        ax.plot(x,z,c='k', lw=2)
        x,y,z=TE[part].T
        ax.plot(x,z,c='k', lw=2)
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('z')
    plt.show()

    #########################################
    #TODO: Get this part workign again

    # # Sorting data points
    # fig2 = plt.figure()
    # ax = fig2.add_subplot(111)

    # for part in parts:
        # for i in range(len(eta_sampled)):
            # eta = eta_sampled[i]
            # print('eta: %.3f, twist: %3f' % (eta_sampled[i], np.degrees(theta_list[i])))
            # rotated_data = am.rotate({'x':data_sampled[part][eta]['x'],
                                      # 'y':data_sampled[part][eta]['z']}, 
                                      # {'x':[],'y':[]},
                                      # {'x':LE_list['x'][i], 'y':LE_list['z'][i]}, 
                                      # theta_list[i], unit_theta = 'rad',
                                      # move_to_origin = True, chord=chord_list[i])
            # ax.scatter(rotated_data['x'],
                       # rotated_data['y'], label=r'$\eta$=%.3f' % eta)        
    # plt.legend()
    # plt.show()
    #pickle.dump(data_sorted, open( "calibrated_values.p", "wb" ) )
