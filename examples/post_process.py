# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 10:36:46 2014

@author: Pedro
"""

from abaqus import *
from abaqusConstants import *
import visualization
from viewerModules import *

import numpy as np

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_displacement(odbName, StepName, components = [2]):

    """
    This ODB reading script does the following:
    -Retrieves the displacement at TIPNODE
    """

    # Open the output database.
    print 'Made it in'
#    odbName = ModelName+'.odb'
#    print odbName
    odb = visualization.openOdb(odbName)
    lastFrame = odb.steps[StepName].frames[-1]
    print 'odb open'


    # Selecting the node(s) to be queried
    pTip = odb.rootAssembly.nodeSets['TIPNODE']

    # Retrieve Y-displacements at the splines/connectors
    print 'Retrieving ALL final displacements at ALL points'
    dispField = lastFrame.fieldOutputs['U']

    print 'Retrieving ALL displacements at TIPNODE'
    dFieldpTip = dispField.getSubset(region=pTip)

    print 'Retrieving only U2 at TIPNODE'
    #Note, U1=data[0], U2=data[1], U3=data[2]
    disppTip = []
    for component in components:
        disppTip.append(dFieldpTip.values[0].data[component-1])
    odb.close()
    if len(disppTip) == 1:
        return disppTip[0]
    else:
        return disppTip


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_coordinates_linear(ModelName,StepName,InstanceName,SetName):

    """
    This ODB reading script does the following:
    -Retrieves the coordinates at SetName for a LINEAR ANALYSIS
	
	DO NOT USE THIS WHEN NLGEOM IS ON! IT WILL DOUBLE DIP!
    """
    Coordinates={'x':[],'y':[]}
    # Open the output database.
    print 'Made it in'
    odbName = ModelName+'.odb'
#    print odbName
    odb = visualization.openOdb(odbName)
    lastFrame = odb.steps[StepName].frames[-1]
    print 'odb open'

    # Selecting the node(s) to be queried
    coordset = odb.rootAssembly.instances[InstanceName.upper()].nodeSets[SetName.upper()]

    # Retrieve Y-displacements at the splines/connectors
    print 'Retrieving ALL final displacements at ALL points'
    dispField = lastFrame.fieldOutputs['U']

    print 'Retrieving ALL displacements at TIPNODE'
    dFieldpTip = dispField.getSubset(region=coordset)

    print 'Retrieving ALL final displacements at ALL points'
    coordField = lastFrame.fieldOutputs['COORD']

    print 'Retrieving ALL displacements at Set'
    cFieldpTip = coordField.getSubset(region=coordset)

    print 'Retrieving coordinates at sets at Set'
    #Note, U1=data[0], U2=data[1], U3=data[2]
    for i in range(len(dFieldpTip.values)):
        Coordinates['x'].append(cFieldpTip.values[i].data[0]+dFieldpTip.values[i].data[0])
        Coordinates['y'].append(cFieldpTip.values[i].data[1]+dFieldpTip.values[i].data[1])
    odb.close()
    return Coordinates

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def findOutputSet(ModelName,StepName,InstanceName,SetName, Output='COORD'):
    output=[]
    time = []
    odbName = ModelName+'.odb'

    odb = visualization.openOdb(odbName)
    for i in range(len(odb.steps[StepName].frames)):
        output_frame = []
        lastFrame = odb.steps[StepName].frames[i]
        print(i, lastFrame.frameValue)
        time.append(lastFrame.frameValue)
        coordset = odb.rootAssembly.instances[InstanceName.upper()].nodeSets[SetName.upper()]

        # Retrieve Y-displacements at the splines/connectors
        dispField = lastFrame.fieldOutputs[Output]

        dFieldpTip = dispField.getSubset(region=coordset)

        for i in range(len(dFieldpTip.values)):
            output_frame.append(dFieldpTip.values[i].data)
        output.append(output_frame)
    odb.close()

    return np.array(output), np.array(time)
    
if __name__ == '__main__':
    import pickle as p
    ModelName = 'fuselage_part_03-new'
    Steps = ['Step-3'] # 'Step-1', 'Step-2', 
    InstanceName = 'Part-3-1'
    SetName = 'Set-1'
    outputNames = ['NT11', 'COORD', 'E', 'SDV2']
    coordinates = {}
    temperatures = {}
    output = {'Time':{}}
    for outputName in outputNames:
        output[outputName] = {}
        for StepName in Steps:
            output_i, time = findOutputSet(ModelName,StepName,InstanceName,
                                         SetName, outputName)
            output[outputName][StepName] = output_i
            output['Time'][StepName] = time                                             
    f = open('outputs.p', 'wb')
    p.dump(output, f)
        