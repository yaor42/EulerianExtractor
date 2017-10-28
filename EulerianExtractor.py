""" 
This script post-processes an Abaqus Odb file and extracts the time 
history of a field output of a body at some fixed spattial location, 
or in mechanics people's word, the Eulerian result of body at one 
spatial location.
 
author: Yao Ren
website: http://yaor.me
license: MIT 
 
Please feel free to adapt this script to your need. But the author would
really appreciate if you keep the above information. Thanks. 
"""
# Define the odb file name, please include .odb extension at the end.
# Type: string
odbName = '2RR_T1-1000_T2-500_N20.odb'

# Define (x, y, z) coordinates of a spatial location of interest. The
# result of interest of body (defined below) at this spatical location
# will be extracted. 
# Type: tupe of floats
poi = (0., 0., 0.0)

# Define the field output name of interest, should not include component.
# Type: string
fieldVarName = 'V'
# Define the component of the vector or tensor type field ouptut, 
# Type: int or SymbolicConstant
# int 1, 2, 3 for vector for now
fieldVarComponent = 1

# Define the instance of model whose result is of interest
# Type: string
instanceName = 'WEB-1'

# Define a potential search range for nodes of the above instance. The poi 
# location should be enclosed in this potential range all the time. It's
# recommonded to define to reduce processing time if model is large. Leave 
# it blank if all nodes of the instance are of interest. Format follows 
# Abaqus nodal path. 
# Type: string 
potentialNodeRange = '4403:6603'        

# End of input

import numpy as np
from random import randint
from odbAccess import *
import visualization
from abaqusConstants import *

odb = openOdb(path=odbName)
instance = odb.rootAssembly.instances[instanceName]
stepRepo = odb.steps

def nLabelStr2IntTuple(nLabelStr):
    """converts an legal nodeList string to a to tuple consists of integer labels of those nodes
    format of nodeList is same as the one used to create nodal path in Abaqus visualization
    example: '1,3:7:2,10:8:-1' will be converted to (1,3,5,7,10,9,8)
    """
    resList = []
    nLists = [nSeg.strip() for nSeg in nLabelStr.split(',')]
    for nSeg in nLists:
        if ':' in nSeg:
            if nSeg.count(':') == 1:
                start, end = [int(x) for x in nSeg.split(':')]
                nTempList = range(start, end+1)
            else:
                start, end, inc = [int(x) for x in nSeg.split(':')]
                if inc < 0:
                    end -= 1
                nTempList = range(start, end, inc)
        else:
            nTempList = [int(nSeg)]
        resList += nTempList
    return tuple(resList)

def createNodeSet(nodeLabels, nodeSetName, instance):
    """create an Abaqus nodeset object, if the nodeset is not created previously
    nodeLabels: a tuple consists of labels of nodes
                or a string
    nodeSetName: a string of nodeSet name, it will be used as the key of repository 
    """
    if nodeSetName not in instance.nodeSets.keys():
        # empty nodeLabesl creates nodeSet of all nodes of a instance
        if nodeLabels == '':
            nodeSet = instance.NodeSet(name=nodeSetName, nodes=instance.nodes[:])
        else:
            if isinstance(nodeLabels, str):
                nodeLabelsTuple = nLabelStr2IntTuple(nodeLabels)
            else:
                nodeLabelsTuple = nodeLabels
            nodeSet = instance.NodeSetFromNodeLabels(name=nodeSetName, nodeLabels=nodeLabelsTuple)
    else:
        nodeSet = instance.nodeSets[nodeSetName]
    return nodeSet        

def getNLabelsCoordSeq(frame, nodeSet):
    nLabelsSeq = ()     # in ascending sequence of node label number
    coordSeq = ()       # in cooresponded sequence of nLabelsSeq
    fieldValues = frame.fieldOutputs['U'].getSubset(region=nodeSet).values
    for v, n in zip(fieldValues, nodeSet.nodes):
        nLabelsSeq += (n.label, )
        coordSeq += (v.data + n.coordinates, )
    
    return nLabelsSeq, coordSeq

def createCoordField(nodeSet, instance):
    lastStep = stepRepo[odb.steps.keys()[-1]]
    if 'COORD' in lastStep.frames[-1].fieldOutputs.keys():
        return
    scratchOdb = session.ScratchOdb(odb=odb)
    for stepKey in stepRepo.keys():
        step = stepRepo[stepKey]
        try:
            ghostStep = scratchOdb.Step(name='Ghost'+stepKey, 
                                        description='Ghost step for your step',
                                        domain=TIME, timePeriod=step.timePeriod)
        except OdbError:
            print 'Ghost step has been created before. Old ghost step will be used.'
            ghostStep = scratchOdb.steps['Ghost'+stepKey]
        if len(ghostStep.frames) != len(step.frames):
            for frame in step.frames:
                ghostFrame = ghostStep.Frame(frameId=frame.frameId, frameValue=frame.frameValue,
                                description='Ghost frame')
                ghostField = ghostFrame.FieldOutput(name='COORD', 
                                       description='initial position plus displacement to get coordinates',
                                       type=VECTOR)
                nLabelsSeq, coordSeq = getNLabelsCoordSeq(frame, nodeSet)
                ghostField.addData(position=NODAL, instance=instance, 
                                   labels=nLabelsSeq , data=coordSeq)
                print ''.join(('Prepare coordinates for ', stepKey, ': ', 
                               str(frame.frameId), '/', str(len(step.frames))))
    return

def closestNodeLabel(nodeSet, poi, step, frame, instance):
    """find the node in nodalRange that has the shortest distance to poi
    return: the int label of the above node
    """
    dist = float('inf')
    if 'COORD' in frame.fieldOutputs.keys():
        coordFieldPotential = frame.fieldOutputs['COORD'].getSubset(region=nodeSet)
    else:
        ghostOdbName = [k for k in session.scratchOdbs.keys() if odbName in k][0]
        coordFieldPotential = session.scratchOdbs[ghostOdbName].steps['Ghost'+step.name].getFrame(frameValue=frame.frameValue).fieldOutputs['COORD']
    poiA = np.array(poi)
    for nCoordValue in coordFieldPotential.values:
        nCoordData = nCoordValue.data
        distNew = np.linalg.norm(poiA-nCoordData)
        if distNew < dist:
            dist = distNew
            closestNLabel = nCoordValue.nodeLabel    
    return closestNLabel

def createXYDataObj(xySequence, xyDataName):
    if xyDataName not in session.xyDataObjects.keys():
        resXYData = session.XYData(data=xySequence, name=xyDataName)
    else:
        resXYData = session.xyDataObjects[xyDataName].setValues(data=xySequence)
    return resXYData

def getVarValue(fieldVarName, fieldVarComponent, frame, regionSet):
    varType = frame.fieldOutputs[fieldVarName].values[-1].type
    if varType == SCALAR:
        pass
    elif varType == VECTOR:
        varValue = frame.fieldOutputs[fieldVarName].getSubset(region=regionSet).values[0].data[fieldVarComponent-1]
    elif varType == TENSOR_2D_PLANAR:
        pass
    elif varType == TENSOR_2D_SURFACE:
        pass
    elif varType == TENSOR_3D_FULL:
        pass
    elif varType == TENSOR_3D_PLANAR:
        pass
    elif varType == TENSOR_3D_SURFACE:
        pass
    return varValue

nodeSetName = 'potentialNodeSetForSearch'
if potentialNodeRange == '' or ' ':
    nodeSetName = 'allInstanceNodes'

nodeSetPotential = createNodeSet(potentialNodeRange, nodeSetName, instance)
createCoordField(nodeSetPotential, instance)

xySeq = ()
if fieldVarName not in stepRepo[stepRepo.keys()[-1]].frames[-1].fieldOutputs.keys():
    print ''.join((fieldVarName, ' is not an Abaqus field output. Please check spelling.'))
else:
    for stepKey in stepRepo.keys():
        step = stepRepo[stepKey]
        frameNum = 0
        stepFrameNum = len(step.frames)
        for frame in step.frames:
            dataPoints = []
            time = step.totalTime + frame.frameValue
            dataPoints.append(time)
            nLabel = closestNodeLabel(nodeSetPotential, poi, step, frame, instance)
            nodeSetName = str(nLabel) + stepKey + str(frame.frameId)
            nodeSet = createNodeSet((nLabel,), nodeSetName, instance)
            #output = frame.fieldOutputs[fieldVarName].getSubset(region=nodeSet).values[0].data[fieldVarComponent-1]
            output = getVarValue(fieldVarName, fieldVarComponent, frame, nodeSet)
            dataPoints.append(output)
            xySeq += (tuple(dataPoints),)
            frameNum += 1
            print ''.join((stepKey, ': ', str(frameNum), '/', str(stepFrameNum)))
    
    spatialXYData = createXYDataObj(xySequence=xySeq, xyDataName='SpatialXYData-'+fieldVarName+str(fieldVarComponent))

