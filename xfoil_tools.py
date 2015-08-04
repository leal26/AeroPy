"""
Created on Sun Mar  9 14:58:25 2014
Last update Fr Jul 20 16:26:40 2015
    
@author: Pedro Leal
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Import necessary modules
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Module that subsistutes the os module and can be used to access the
# command line 
import subprocess as sp

# To check for already existing files and delete them
import os

# Numpy module is imported in case one of the inputs is a numpy array 
# and for mathematical analysis
import numpy as np

# Modules necessary for saving multiple plots
import shutil

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       	Core Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def call(airfoil, alfas='none', output='Cp', Reynolds=0, Mach=0, plots=False,
         NACA=False, GDES=False,iteration=10):
    """ Call is a function that enables the use of xfoil through 
        Python.
    The input variables are:
        - airfoil: if NACA is false, airfoil is the name of the plain 
          filewhere the airfoil geometry is stored (variable airfoil).
          If NACA is True, airfoil is the naca series of the airfoil 
          (i.e.: naca2244). By default NACA is False.
          
        - alfas: list/array/float/int of angles of attack.
        
        - output: defines the kind of output desired from xfoil. There 
          are four posssible choices:
              - Cp: generates files with Pressure coefficients for
                desired alfas.
              - Dump: generates file with Velocity along surface, Delta
                star,theta and Cf vs s,x,y for several alfas.
              - Polar: generates file with CL, CD, CM, CDp, Top_Xtr, 
                Bot_Xtr.
		      - Alfa_L_0: generates a file with the value of the angle of
			    attack that lift is equal to zero.
              - Coordinates: returns the coordinates of a NACA airfoil.
              
          By default, Cp is chosen.
              
        - Reynolds: Reynolds number in case the simulation is for a 
          viscous flow. In case not informed, the code will assume
          inviscid.
          
        - Mach: Mach number in case the simulation has to take in 
          account compressibility effects through the Prandtl-Glauert
          correlation. If not informed, the code will not use the
          correction. For logical reasons, if Mach is informed a 
          Reynolds number different from zero must also be informed.
          
        - plots: the code is able to save in a .ps file all the plots
          of Cp vs.alfa. By default, this option is deactivated.
          
        - NACA: Boolean variable that defines if the code imports an 
          airfoil from a file or generates a NACA airfoil.
          
        - GDES: XFOIL function that improves the airfoil shape in case
          the selected points do not provide a good shape. The CADD
          function is also used. For more information about these
          functions, use the XFOIL manual.
          
        - iteration: changes how many times XFOIL will try to make the
          results converge. Speciallt important for viscous flows
          
    As a side note, it is much more eficient to run a single run with
    multiple angles of attack rather than multiple runs, each with a 
    single angle of attack

    Created on Sun Mar  9 14:58:25 2014
    Last update Fr Jul 13 15:38:40 2015
    
    @author: Pedro Leal (Based on Hakan Tiftikci's code)
    """
    

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                               Functions
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
    def issueCmd(cmd,echo=True):
        """Function that submits commands through the PIPE to the 
        command line, therefore leading the commands to xfoil.
        """
        
        ps.stdin.write(cmd+'\n')
        if echo:
            print cmd

            
    def submit(output,alfa):
        """Submits job to xfoil and saves file.
        
        Standard output file= function_airfoil_alfa.txt, where alfa has
        4 digits, where two of them are for decimals. i.e. 
        cp_naca2244_0200. Analysis for Pressure Coefficients for a 
        naca2244 at an angle of degrees.
        
        Possible to output other results such as theta, delta star 
        through the choice of the ouput, but not implemented here.
        """
        
        if output=="Alfa_L_0":
			issueCmd('CL 0')
			
        else:
		# Submit job for given angle of attack
		issueCmd('ALFA %.4f' % (alfa,))
		if plots==True:
			issueCmd('HARD')
			shutil.copyfile('plot.ps','plot_{!s}_{!s}_{!s}.ps'.format(output,
							airfoil,alfa))
		if output=='Cp':    
			# Creating the file with the Pressure Coefficients
			filename=file_name(airfoil,alfas,output)
			try:
				os.remove(filename)
			except OSError:
				pass 
			#issueCmd('CPWR {!s}'.format(filename))
			issueCmd('CPWR %s' % filename)

		if output=='Dump':    
			# Creating the file with the Pressure Coefficients
			filename=file_name(airfoil,alfas,output)
			try:
				os.remove(filename)
			except OSError:
				pass      
	
			issueCmd('DUMP %r' % filename)
     
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                Characteristics of the simulation
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # By default the code considers the flow to be inviscid.
    Viscid=False    
    if Reynolds!= 0:
        Viscid=True
    # Is alpha given or not?(in case of Alfa_L_0, then alfas=False)
    if alfas!='none':
        print type(alfas)
        # Single or multiple runs?
        if type(alfas)==list or type(alfas)==np.ndarray:
            Multiple=True
        elif type(alfas)==int or type(alfas)==float:
            Multiple=False
    elif output=="Alfa_L_0" and alfas=='none':
        Multiple=False
    elif output=="Alfa_L_0" and alfas!='none':
        raise Exception("To find alpha_L_0, alfas must not be defined")
    elif output!="Alfa_L_0" and alfas=='none':
        raise Exception("To find anything except alpha_L_0, you need to "
		    "define the values for alfa")	
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                           Start Xfoil
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    """For communication with the xfoil through the command line the 
    Popen class from subprocess is used. stdin and stdout both 
    represent inputs and outputs with the process run on the command 
    line, in this case, xfoil.
    
    class Popen(args, bufsize=0, executable=None,
                stdin=None, stdout=None, stderr=None,
                preexec_fn=None, close_fds=False, shell=False,
                cwd=None, env=None, universal_newlines=False,
                startupinfo=None, creationflags=0):
    """
    # Calling xfoil with Poper 
    ps = sp.Popen(['xfoil.exe'],
                  stdin=sp.PIPE,
                  stdout=None,
                  stderr=None)
                  
    # Loading geometry
    if NACA==False:
        #issueCmd('load {!s}'.format(airfoil))
        issueCmd('load %s' % airfoil)
    else:
        #issueCmd('{!s}'.format(airfoil))
        issueCmd('%s' % airfoil)

    # If output equals Coordinates, no analysis will be realized, only the 
    # coordinates of the shape will be outputed        
    if output=='Coordinates':
        issueCmd('SAVE')
        issueCmd(airfoil + '.txt')
        # In case there is alread a file with that name, it will replace it.
        # The yes stands for YES otherwise Xfoil will do nothing with it.
        issueCmd('Y')
    else:
        # Once you load a set of points in Xfoil you need to create a 
        # name, however we do not need to give it a name
        issueCmd('')
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #             Adapting points for better plots
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         
        """
        The following part of the code was not implemented in hope that
        the Cp be given for the given airfoil points, not for points 
        created by xfoil.    
        
        GDES   (enter GDES menu)        |
        CADD   (add points at corners)  |  These commands are optional,
               (accept default input)   |  and are recommended only for
               (accept default input)   |  Eppler and Selig airfoils
               (accept default input)   |  to give smoother LE shapes
               (return to Top Level)    |
         
        PANEL  (regenerate paneling since better panel node spacing is
                needed)
        """
        if GDES==True:
            issueCmd('GDES')  # enter GDES menu
            issueCmd('CADD')  # add points at corners
            issueCmd('')      # accept default input
            issueCmd('')      # accept default input
            issueCmd('')      # accept default input
            issueCmd('')      # accept default input
            issueCmd('PANEL') # regenerate paneling
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Opening OPER module in Xfoil
        issueCmd('OPER')
    
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #                Applying effects of vicosity
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        issueCmd('iter')
        issueCmd('%d' % iteration)
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #                Applying effects of vicosity
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
        if Viscid==True:
            # Defining the system as viscous
            issueCmd('v')
            # Defining Reynolds number
            issueCmd('%f' % Reynolds)
    
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #    Defining Mach number for Prandtl-Gauber correlation
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        #issueCmd('MACH {!s}'.format(Mach))
        issueCmd('MACH %s' % Mach)
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #                     Submitting
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if output=='Polar' or output=='Alfa_L_0':
            issueCmd('PACC')
            # All file names in this library are generated by the
            # filename functon.
            filename=file_name(airfoil,alfas,output)
            try:
                os.remove(filename)
            except OSError:
                pass
                
            issueCmd('%s' % filename)
            issueCmd('')
        
        # For several angles of attack
        if Multiple==True:
            for alfa in alfas:
                submit(output,alfa)
            
        # For only one angle of attack
        if Multiple==False:
            submit(output,alfas)
            
        # Exiting
        # From OPER mode
        issueCmd('')
    # From xfoil
    issueCmd('QUIT')
    # From stdin
    ps.stdin.close()
    # From popen
    ps.wait()
    
def create_x(c):
    """
    The ojective of this function is to create the set of points along
    the chord. The x output is conviniently ordered from TE to LE. 
    
    Because Xfoil it is not efficient or sometimes possible to create a
    uniform distribution of points because the front part of the 
    airfoil requires a big amount of points to create a smooth surface 
    representing the round leading edge. However there is a maximum of 
    points that Xfoil will accept, therefore there is a need to 
    overcome this obstacle. This was done by dividing the airfoil in 
    to 3 parts.
    
        - tip: correpondent to 0 to .3% of the chord length, it is the
          most densily populated part of the airfoil to compensate the
          wide variation of slopes
        
        - middle: correpondent to .3% to 30% of the chord length (
          Shortly above a quarter of the chord length). The second most
          densily populated and the segment with the most amount of 
          points. Represents the round section of the airfoil except 
          for the tip. Such an amount of points is necessary because it
          is where most of the lift is generated.
        
        - endbody: correspondent to 30% to 100% of the chord length. 
          The less densily populated section. Such unrefined mesh is
          possible because of the simplistic geometry of endbody's
          airfoil, just straight lines.
        
    Created on Thu Feb 27 2014
    
    @author: Pedro Leal
    """
    
    max_point=c/4
    limit=max_point+0.05*c    
    nose_tip=0.003*c

    # Amount of points for each part
    N_tip=10
    N_middle=180
    N_endbody=40

    x_endbody=np.linspace(c,limit,N_endbody)  
    x_middle=np.linspace(limit,nose_tip,N_middle)
    x_tip=np.linspace(nose_tip,0,N_tip)
    
    # Organizing the x lists in a unique list without repeating any
    # numbers
    x2=np.append(x_middle,np.delete(x_tip,0))
    x=np.append(x_endbody,np.delete(x2,0))
    return x


def create_input(x,y_u,y_l,filename):
    """
    Function generates a plain file that xfoil can read.
    
    XFOIL only reads file from the TE to the LE from the upper part
    first and then from the LE to the TE through the pressure surface.
    
    Inputs:
        - x: list of coordinates along the chord
        
        - y_u: list of coordinates normal to the chord for the upper 
          surface
        
        - y_l: list of coordinates normal to the chord for the lower
          surface
        
        - file_name: label used for the file created
        
    Created on Thu Feb 27 2014
    
    @author: Pedro Leal"""
    
    # XFOIL likes to read the files from the TE to the LE from the
    # upper part first and then from the LE to the TE through the
    # pressure surface
    x_upper=x
    x_under=np.delete(x_upper,-1)[::-1]
    x=np.append(x_upper,x_under)
    
    y_l=np.delete(y_l,-1)[::-1]
    y=np.append(y_u,y_l)
    
    # Creating files for xfoil processing
    DataFile = open(filename,'w')
    for i in range(0,len(x)):
        DataFile.write('     %f    %f\n' % (x[i],y[i]))
    DataFile.close()
    return 0
    
def prepare_xfoil(Coordinates_Upper, Coordinates_Lower, chord, 
                  reposition=False, FSI=False):
    # The upper and lower functions will be the points in ordered 
    # fashion. Because of the way that XFOIL works the points start at
    # the Trailing Edge on the upper surface going trough the Leading 
    # Edge and returning to the Trailing Edge form the bottom surface.
    
    def Reposition(CoordinatesU,CoordinatesL):

        # Find the coordinates of the trailing edge
        LE={'x':0,'y':0}
        cx=CoordinatesU['x']
        LE['x']=min(cx)
        index_LE=cx.index(LE['x'])
        LE['y']=cx[index_LE]
        
        All_Rotated_Coordinates={}
        count=0
        for Coordinates in [CoordinatesU,CoordinatesL]: 
            # Repositioning Leading Edge
            for key in Coordinates:
                c=Coordinates[key]
                c=[i - LE[key] for i in c]
                Coordinates[key]=c
                
        # Find the coordinates of the trailing edge. Because of the 
        # thickness of the TE, it is necessary to find the point with
        # max(x) for both surfaces and take the average between them
        # to find the actual TE.     
        TE={'x':0,'y':0}
        
        cxU=CoordinatesU['x']
        cyU=CoordinatesU['y']
        TExU=max(cxU)
        index_TE=cxU.index(TExU)
        TEyU=cyU[index_TE]
        
        cxL=CoordinatesL['x']
        cyL=CoordinatesL['y']
        TExL=max(cxL)
        index_TE=cxL.index(TExL)
        TEyL=cyL[index_TE]
        
        TE['x']=(TExU+TExL)/2.
        TE['y']=(TEyU+TEyL)/2.
        
        # Rotating according to the position of the trailing edge
        theta=math.atan(TE['y']/TE['x'])
        
        # Rotation transformation Matrix
        T=[[math.cos(theta),math.sin(theta)],
            [-math.sin(theta),math.cos(theta)]]
            
        for Coordinates in [CoordinatesU,CoordinatesL]:         
            Rotated_Coordinates={'x':[],'y':[]}
            for i in range(len(Coordinates['x'])):
                cx=Coordinates['x'][i]
                cy=Coordinates['y'][i]
                rot_x=T[0][0]*cx+T[0][1]*cy
                rot_y=T[1][0]*cx+T[1][1]*cy
                Rotated_Coordinates['x'].append(rot_x)
                Rotated_Coordinates['y'].append(rot_y)
            All_Rotated_Coordinates['%s' % count]=Rotated_Coordinates
            count+=1
        return All_Rotated_Coordinates['0'],All_Rotated_Coordinates['1']
        
    upper=[]
    lower=[]
    print "Starting to prepare points"
        
    # At first we'll organize the files by its x values
    for i in range(len(Coordinates_Upper['x'])):
        # For each x value, we will check the correpondent y value so
        # that we can classify them as upper or lower
        upper.append([Coordinates_Upper['x'][i]/chord,
                      Coordinates_Upper['y'][i]/chord])

    for i in range(len(Coordinates_Lower['x'])):
        # For each x value, we will check the correpondent y value so 
        # that we  can classify them as upper or lower
        lower.append([Coordinates_Lower['x'][i]/chord,
                      Coordinates_Lower['y'][i]/chord])
    print "Sorting Stuff up"

    if reposition == True:
        # Sort in a convenient way for calculating the error
        upper=sorted(upper,key=lambda coord:coord[0], reverse=False)
        lower=sorted(lower,key=lambda coord:coord[0], reverse=False)
        print 'Repositioning'
        cu = {'x':[], 'y':[]}
        cl = {'x':[], 'y':[]}
        for i in range(len(upper)):
            cu['x'].append(upper[i][0])
            cu['y'].append(upper[i][1])
        for i in range(len(lower)):
            cl['x'].append(lower[i][0])
            cl['y'].append(lower[i][1])
        upper, lower = Reposition(cu,cl)
        print "Done preparing points"
        return upper,lower
    elif FSI == True:
        upper = sorted(upper, key=lambda coord:coord[0], reverse=False)
        lower = sorted(lower, key=lambda coord:coord[0], reverse=False)   
        print "Done preparing points"
        return upper,lower
    else:
        # Sort in a way that comprehensible for xfoil and elimates the
        # repeated point at the LE
        upper = sorted(upper, key=lambda coord:coord[0], reverse=True)
        lower = sorted(lower, key=lambda coord:coord[0], reverse=False)[1:]
        Coordinates = upper + lower
        print "Done preparing points"
        return Coordinates
    
def output_reader(filename, separator='\t', output=None, rows_to_skip=0, header=0):
    """Function that opens files of any kind. Able to skip rows and
    read headers if necessary.
    
    Inputs:
        - filename: just the name of the file to read.
        
        - separator: Main kind of separator in file. The code will
          replace any variants of this separator for processing. Extra
          components such as \n kg m are all eliminated.
          
        - output: defines what the kind of file we are opening to
          ensure we can skip the right amount of lines. By default it 
          is None so it can open any other file.
          
        - rows_to_skip: amount of rows to initialy skip in the file. If
          the output is different then None, for the different types of
          files it is defined as:
          - Polar files = 10
          - Dump files = 0
          - Cp files = 2
          
        - header: The header list will act as the keys of the output
          dictionary. For the function to work, a header IS necessary.
          If not specified by the user, the function will assume that
          the header can be found in the file that it is opening.
    Output:
    Dictionary with all the header values as keys
    
    Created on Thu Mar 14 2014
    
    @author: Pedro Leal
    """
    # In case we are using an XFOIL file, we define the number of rows
    # skipped
    if output=='Polar' or output=='Alfa_L_0':
        rows_to_skip=10
    elif output=='Dump':
        rows_to_skip=0
    elif output=='Cp':
        rows_to_skip=2
    # n is the amount of lines to skip
    Data={}
    header_done=False
    count_skip=0
    
    with open (filename, "r") as myfile:
        # Jump first lines which are useless
        for line in myfile:
            if count_skip<rows_to_skip:
                count_skip+=1
                # Basically do nothing
            elif header_done==False:
                # If the user did not specify the header the code will
                # read the first line after the skipped rows as the
                # header
                if header==0:
                    # Open line and replace anything we do not want (
                    # variants of the separator and units)
                    line=line.replace(separator+separator+separator+separator+
                    separator+separator,' ').replace(separator+separator+
                    separator+separator+separator,' ').replace(separator+
                    separator+separator+separator,' ').replace(separator+
                    separator+separator,' ').replace(separator+separator,
                    ' ').replace("\n","").replace("(kg)","").replace("(m)",
                    "").replace("(Pa)","").replace("(in)","").replace("#","")
                    
                    header=line.split(' ')
                    n_del=header.count('')
                    for n in range(0,n_del):
                        header.remove('')
                    for head in header:
                        Data[head]=[]
                    # To avoid having two headers, we assign the False 
                    # value to header which will not let it happen
                    header_done=True
                # If the user defines a list of values for the header, 
                # the code reads it and creates libraries with it.
                elif type(header)==list:
                    for head in header:
                        Data[head]=[]
                    header_done=True
            else:
                line=line.replace(separator+separator+separator,
                ' ').replace(separator+separator,' ').replace(separator,' ').replace("\n",
                "").replace('---------','').replace('--------',
                '').replace('-------','').replace('------','').replace('-',
                ' -')
                
                line_components=line.split(' ')    
                
                n_del=line_components.count('')
                for n in range(0,n_del):
                    line_components.remove('')
                   
                if line_components!=[]:
                    for j in range(0,len(line_components)):
                        Data[header[j]].append(float(line_components[j]))
                # else DO NOTHING!
    return Data
    
def alfa_for_file(alfa):
    alfa='%.2f' % alfa
    inter,dec=alfa.split('.')
    inter_number=int(inter)
    inter='%.2d' % inter_number
    if inter_number<0:
        inter='n'+inter
    alfa=inter+dec
    return alfa
    
def file_name(airfoil,alfas='none',output='Cp'):
    """Function creates name of files generated by XFOIL.
    
    The input variables are:
        - airfoil: the name of the plain file where the airfoil 
          geometry is stored (variable airfoil).
          
        - alfas: list/array/float/int of a single angle of attack for
          Cp and Dump, but the whole list for a Polar. Only the initial
          and the final values are used
        
        - output: defines the kind of output desired from xfoil. There
          are three posssible choices:
              - Cp: generates files with Pressure coefficients for 
                desired alfas
              - Dump: generates file with Velocity along surface, Delta
                star and theta and Cf vs s,x,y for several alfas
              - Polar: generates file with CL, CD, CM, CDp, Top_Xtr,
                Bot_Xtr
			  - Alfa_L_0: calculate the angle of attack that lift is
                zero
                
    The output has the following format:
        - for Cp and Dump: output_airfoil_alfa
           i.e. Cp_naca2244_0200 (alfa =2.00 degress)
        - for Polar: Polar_airfoil_alfa_i_alfa_f
           i.e. Polar_naca2244_n0200_0200 (alfa_i=-2.00 degrees and 
                                           alfa_f=2.00)
		- for Alpha_L_0: Alpha_L_0_airfoil
		   i.e. Alpha_L_0_naca2244
            
    By default, Cp is chosen.
            
    Created on Thu Mar 16 2014
        
    @author: Pedro Leal 
    """
	# At first verify if alfas was defined
    if alfas=='none':
		filename='%s_%s' % (output,airfoil)
    elif alfas!='none':
		if output=='Cp' or output=='Dump':
			if type(alfas)==list:
				alfas=alfas[0]
			alfa=alfa_for_file(alfas))

			filename='%s_%s_%s' % (output,airfoil,alfa)
				
		if output=='Polar':
			# In case it is only for one angle of attack, the same
            # angle will be repeated. This is done to keep the 
            # formating
			if type(alfas)==int or type(alfas)==float:
				alfas=[alfas]
			alfa_i=alfa_for_file(alfas[0])
			alfa_f=alfa_for_file(alfas[-1])
			# Name of file with polar information

			filename='%s_%s_%s_%s' % (output,airfoil,alfa_i,alfa_f)
    return filename
    
def M_crit(airfoil,pho,speed_sound,lift,c):
    M_list=np.linspace(0.3,0.7,20)
    alfas=np.linspace(-15,5,21)
    Data_crit={}
    Data_crit['M']=0
    Data_crit['CL']=0
    Data_crit['alpha']=0
    for M in M_list:
        cl=(np.sqrt(1-M**2)/(M**2))*2*lift/pho/(speed_sound)**2/c
        call(airfoil,alfas,output='Polar',NACA=True)
        filename=file_name(airfoil,alfas,output='Polar')
        Data=output_reader(filename,' ',10)
        previous_iteration=Data_crit['CL']
        for i in range(0,len(Data['CL'])):
            if Data['CL'][i]>=cl and M>Data_crit['M']:
                print M
                Data_crit['M']=M
                Data_crit['CL']=Data['CL'][i]
                Data_crit['alpha']=Data['alpha'][i]
#        if Data_crit['CL']==previous_iteration:
    return Data_crit

def air_properties(height,unit='feet'):
    """ Function to calculate air properties for a given height (m or ft).
    
    Sources: 
      - http://en.wikipedia.org/wiki/Density_of_air#Altitude
      - http://aerojet.engr.ucdavis.edu/fluenthelp/html/ug/node337.htm
      
    Created on Thu May 15 14:59:43 2014
    
    @author: Pedro Leal
    """
    # height is in m
    if unit=='feet':
        height=0.3048*height
    elif unit != 'meter':
        raise Exception('air_properties can onlu understand feet and meters')
        
    #==================================================================
    # Constants
    #==================================================================
    # Sea level standard atmospheric pressure
    P0=101325. # Pa
    
    # Sealevel standard atmospheric temperature
    T0=288.15 # K
    
    # Earth-surface gravitational acceleration
    g=8.80655 # m/s2
    
    # Temperatur lapse rate, 0.0065 K/m
    L=0.0065 # K/m
    
    # Ideal (Universal) gas constant
    R=8.31447 # J/(mol K)
    
    # Molar mass of dry air
    M=0.0289644 #kg/mol
    
    # Specific R for air
    R_air=R/M
    # Sutherland's law coefficients
    C1=1.458e-6 #kg/m.s.sqrt(K)
    C2=110.4 #K
    
    #==================================================================
    # Temperature
    #==================================================================
    #Temperature at altitude h meters above sea level is approximated 
    # by the following formula (only valid inside the troposphere):
    
    T = T0 - L*height
    
    #==================================================================
    # Pressure
    #==================================================================
    P=P0*(1-L*height/T0)**(g*M/(R*L))
    
    #==================================================================
    # Density 
    #==================================================================
    density=P*M/(R*T)
    
    #==================================================================
    # Dynamic Viscosity (Sutherland equation with two constants)   
    #==================================================================
    dyn_viscosity=(C1*T**(3./2))/(T+C2)
    
    return {'Density':density,'Dynamic Viscosity':dyn_viscosity,
            'Atmospheric Temperature':T,'R air':R_air, 
            'Atmospheric Pressure': P}

def Reynolds(height,V,c):
    """Simple function to calculate Reynolds for a given height.

    @author: Pedro Leal
    Created in Jul 17 2015    
    """
    Air_Data=air_properties(height,unit='feet')
    rho=Air_Data['Density']
    L=c
    nu=Air_Data['Dynamic Viscosity']
    return rho*V*L/nu        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       	Utility functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def find_coefficients(airfoil,alpha,Reynolds=0,iteration=10,NACA=True):
    filename=file_name(airfoil,alpha,output='Polar')
    # If file already exists, there is no need to recalculate it.
    if not os.path.isfile(filename):
		call(airfoil,alpha,Reynolds=Reynolds,output='Polar',iteration=10,
           NACA=NACA)
    coefficients={}
    # Data from file    
    Data=output_reader(filename,output='Polar')
    for key in Data:
        coefficients[key] = Data[key][0]
    return coefficients
	
def find_alpha_L_0(airfoil, Reynolds=0, iteration=10, NACA=True):
    filename = file_name(airfoil,output='Alfa_L_0')
    # If file already exists, there no need to recalculate it.
    if not os.path.isfile(filename):
        call(airfoil,output='Alfa_L_0',NACA=NACA)
    alpha = output_reader(filename,output='Alfa_L_0')['alpha'][0]
    return alpha
    
def force_shell(Data,chord, half_span, height, Velocity, thickness=0, txt=False):
    # Height is in feet
    # If the Shell is an extrude, it needs to take in consideration
    # that there is a skin thickness outwards of the outer mold.
    # If the Shell is na planar, there is no need for such a 
    # consideration
    Air_properties = air_properties(height,unit='feet')
    atm_pressure = Air_properties['Atmospheric Pressure']
    air_density = Air_properties['Density']
    if thickness == 0:
        Data['Force'] = map(lambda Cp:(Cp*0.5*air_density*(Velocity**2)+atm_pressure)*chord*half_span,Data['Cp'])
        Data['x'] = map(lambda x:(chord)*x,Data['x'])
        Data['y'] = map(lambda x:(chord)*x,Data['y'])
    else:
        Data['Force'] = map(lambda Cp:(Cp*0.5*air_density*(Velocity**2)+atm_pressure)*chord*half_span,Data['Cp'])
        Data['x'] = map(lambda x:(chord-2.*thickness)*x+thickness,Data['x'])
        Data['y'] = map(lambda x:(chord-2.*thickness)*x,Data['y'])
    Data['z'] = [0] * len(Data['x'])
    
    PressureDistribution=zip(Data['x'],Data['y'],Data['z'],Data['Force'])
#    elliptical_distribution=np.sqrt(1.-(Data['z']/half_span)**2) 
#    if txt==True:
#        DataFile = open('Force_shell.txt','w')
#        DataFile.close()
#        for j in range(N):
#            for i in range(len(Data['x'])):
#                    DataFile = open('Force_shell.txt','a')
#                    DataFile.write('%f\t%f\t%f\t%f\n' % (
#                        Data['x'][i],
#                        Data['y'][i],
#                        Data['z'][j],
#                        elliptical_distribution[j]*Data['Force'][i]))
#                    DataFile.close()
#        return 0
#    else:
#        PressureDistribution=()
#        for j in range(N):
#            for i in range(len(Data['x'])):
#                    PressureDistribution=PressureDistribution+((Data['x'][i],
#                        Data['y'][i],Data['z'][j],
#                        elliptical_distribution[j]*Data['Pressure'][i]),)
    return PressureDistribution

def pressure_shell(Data,chord, thickness, half_span, air_density, Velocity, N, txt=False,
                   llt_distribution=False, distribution='Elliptical'): 
                       
    Data['Pressure'] = map(lambda Cp:Cp*0.5*air_density*(Velocity**2)*chord,Data['Cp'])
    Data['x'] = map(lambda x:(chord-2.*thickness)*x+thickness,Data['x'])
    Data['y'] = map(lambda x:(chord-2.*thickness)*x,Data['y'])
    DataFile = open('Pressure_shell.txt','w')
    DataFile.close()
    if distribution=='Elliptical':
        Data['z']=np.linspace(0,half_span,N)
        distribution=np.sqrt(1.-(Data['z']/half_span)**2) 
    elif distribution=='LLT':
        Data['z']=np.linspace(0,half_span,len(distribution))
        distribution=llt_distribution
    if txt==True:
        for j in range(N):
            for i in range(len(Data['x'])):
                    DataFile = open('Pressure_shell.txt','a')
                    DataFile.write('%f\t%f\t%f\t%f\n' % (
                        Data['x'][i],
                        Data['y'][i],
                        Data['z'][j],
                        distribution[j]*Data['Pressure'][i]))
                    DataFile.close()
        return 0
    else:
        PressureDistribution=()
        for j in range(N):
            for i in range(len(Data['x'])):
                    PressureDistribution=PressureDistribution+((Data['x'][i],
                        Data['y'][i],Data['z'][j],
                        distribution[j]*Data['Pressure'][i]),)
        return PressureDistribution
        
def pressure_shell_2D(Data,chord,thickness,half_span,height,Velocity,N,txt=False):
    Air_properties=air_properties(height,unit='feet')
    air_density=Air_properties['Density']
    
    Data['Pressure']=map(lambda Cp:Cp*0.5*air_density*(Velocity**2)*chord,Data['Cp'])
    Data['x']=map(lambda x:(chord-2.*thickness)*x+thickness,Data['x'])
    Data['y']=map(lambda x:(chord-2.*thickness)*x,Data['y'])
    DataFile = open('Pressure_shell.txt','w')
    DataFile.close()
    Data['z']=np.linspace(0,half_span,N)
    if txt==True:
        for j in range(N):
            for i in range(len(Data['x'])):
                    DataFile = open('Pressure_shell.txt','a')
                    DataFile.write('%f\t%f\t%f\t%f\n' % (
                        Data['x'][i],
                        Data['y'][i],
                        Data['z'][j],
                        Data['Pressure'][i]))
                    DataFile.close()
        return 0
    else:
        PressureDistribution=()
        for j in range(N):
            for i in range(len(Data['x'])):
                    PressureDistribution=PressureDistribution+((Data['x'][i],
                        Data['y'][i],Data['z'][j],Data['Pressure'][i]),)
        return PressureDistribution