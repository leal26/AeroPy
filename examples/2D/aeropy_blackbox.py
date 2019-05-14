import numpy as np
import array

from aeropy.geometry.airfoil import CST, create_x, intersect_curves
from aeropy.aero_module import Reynolds, air_properties
from aeropy.xfoil_module import find_coefficients, create_input
from aeropy.morphing.camber_2D import calculate_dependent_shape_coefficients, calculate_strains, plot_airfoil

def aerodynamic_performance(AC, psi_spars, Au_P, Al_P, c_P, deltaz, alpha, H, V):
    morphing_direction = 'forwards'
    
    air_data = air_properties(H, unit='feet')
    density = air_data['Density']
    dyn_pressure = .5*density*V**2
    
    # Generate dependent shape coefficients
    # try:
    Au_C, Al_C, c_C, spar_thicknesses = calculate_dependent_shape_coefficients(
                                                        AC,
                                                        psi_spars, Au_P, Al_P,
                                                        deltaz, c_P, morphing=morphing_direction)

    # print 'Reynolds: ', Reynolds(H, V, c_C)
    # Generate aifoil file
    airfoil = 'test'
    x = create_x(1., distribution = 'linear',n=300)
    y = CST(x, 1., [deltaz/2., deltaz/2.], Au = Au_C, Al= Al_C)

    # Get strain data
    strains, av_strain = calculate_strains( Au_P, Al_P, c_P, Au_C, Al_C, c_C, deltaz, psi_spars, spar_thicknesses)
                    
    intersections = intersect_curves(x, y['l'], x, y['u'])
    print(intersections,intersections[0][1:])
    if len(intersections[0][1:]) == 0:
        # print y
        create_input(x, y['u'], y['l'], airfoil, different_x_upper_lower = False)

        # Get aerodynamic data
        print(airfoil, alpha, Reynolds(H, V, c_C))
        Data = find_coefficients(airfoil, alpha, Reynolds=Reynolds(H, V, c_C), 
                                 iteration=200, NACA=False, delete=True, PANE = True, GDES=True)


        # plot_airfoil(AC, psi_spars, c_P, deltaz, Au_P, Al_P, image = 'save', iteration=counter, dir = airfoil+'_dir')
              
        # filtering data (for now I only care about negative strains
        str_output = {'CL':Data['CL'], 'CD':Data['CD'], 'CM':Data['CM'], 
                          'av_strain':av_strain, 'Au_C':Au_C, 'Al_C': Al_C, 'spars':psi_spars}

        if Data['CM']==None:
            str_output['lift'] = None
            str_output['drag'] = None
            str_output['moment'] = None
        else:
            str_output['lift'] = Data['CL']/dyn_pressure/c_C,
            str_output['drag'] = Data['CD']/dyn_pressure/c_C,
            str_output['moment'] = Data['CM']/dyn_pressure/c_C
        for i in range(len(strains)):
            str_output['strain_'+str(i)] = strains[i]

        # Writing to a text file
        # f_worker = open(str(airfoil) + '.txt', 'wb')   
        # for i in range(len(key_list)):
            # if i != len(key_list)-1:
                # if key_list[i][:1] == 'Au':
                    # f_worker.write('%f\t' % str_output[key_list[i][:1]+'_C'][int(key_list[i][-1])]) 
                # else:
            # else:
                # if key_list[i][:1] == 'Au':
                    # f_worker.write('%f\n' % str_output[key_list[i][:1]+'_C'][int(key_list[i][-1])]) 
    else:
        str_output = {'CL':1000, 'CD':None, 'CM':None, 'av_strain':av_strain, 'spars':psi_spars,
                        'Au_C':Au_C, 'Al_C': Al_C, 'lift': None, 'drag': None, 'moment':None}
        for i in range(len(strains)):
            str_output['strain_'+str(i)] = strains[i]
    # except:
        # str_output = {'CL':None, 'CD':None, 'CM':None, 'av_strain':None,
                        # 'Au_C':[None]*len(AC), 'Al_C': [None]*len(AC), 
                        # 'lift': None, 'drag': None, 'moment':None,
                        # 'strains':[None]*(len(AC)-1), 'spars':psi_spars}     
    return  str_output

if __name__ == '__main__':
    #==============================================================================
    # Parameters
    #==============================================================================
    # morphing_direction = 'forwards'

    # # Velocity
    # V = 30

    # # Angle of attack
    # alpha = 0.

    # # Altitude (ft)
    # H = 10000

    # # Parameter
    # c_P = 1.                  #m
    # deltaz = 0.*c_P    #m

    # # NACA0012
    # # Au_P = [0.16718529, 0.14117793,0.14383898]
    # # Au_P = [0.172802, 0.167353, 0.130747, 0.172053, 0.112797, 0.168891]
    # Au_P = array.array('d',[0.1828, 0.1179, 0.2079, 0.0850, 0.1874])
    # Al_P = Au_P

    # #==============================================================================
    # # Training
    # #==============================================================================

    # AC = array.array('d',[0.1,0.1,0.1,0.1])
    # psi_spars = array.array('d',[0.1,0.3,0.6,0.8])
    # print aerodynamic_performance()
    
    f = open('inputs.txt','r')
    line = f.read()
    line = line.split('\t')
    f.close()
	
    inputs = []
    for i in range(len(line)-1):
        inputs.append(float(line[i]))
    
    AC =  [0.370067,0.397402,0.328316,0.167719]
    psi_spars = inputs[4:8]
    Au_P = inputs[8:13]
    Al_P = inputs[13:18]
    c_P = inputs[18]
    deltaz = inputs[19]
    alpha = 7.1111
    H = inputs[21]
    V = 22.6667
    
    print(AC, alpha, H, V )  
    output = aerodynamic_performance(AC, psi_spars, Au_P, Al_P, c_P, deltaz, alpha, H, V)
    f = open('outputs.txt','w')

    print(output)
    if output['CL'] != None and output['CD'] != None and output['CM'] != None:
        f.write('%6.5f\t%6.5f\t%6.5f\t' % (output['CL'], output['CD'], output['CM']))
    else:
        f.write('%6.5f\t%6.5f\t%6.5f\t' % (0, 10, 0))
    f.close()