""" 
This script post-processes an Abaqus Odb file and extracts the time 
history of a field output of a body at some fixed spattial location, 
or in mechanics people's word, the Eulerian result of body at one 
spatial location.

The method used in this script is to find the node of the instance
input that has the closest distance to the spatial location input,
and extract requested field ouptut value from Odb. The process is
repeated for all frames of all steps of the Odb.
 
author: Yao Ren
website: http://yaor.me
license: MIT 
 
Please feel free to adapt this script to your need. But the author would
really appreciate if you keep the above information. Thanks. 
"""
# Define the odb file name, please include .odb extension at the end.
# Type: string
odbName = '2RR_T1-1000_T2-500_N20-FullIntg.odb'

# Define (x, y, z) coordinates of a spatial location of interest. The
# result of interest of body (defined below) at this spatical location
# will be extracted. 
# Type: tupe of floats
poi = (-1, 0, 0)

# Define the field output name of interest, should not include component.
# Type: string
fieldVarName = 'S'

# Define the component of the vector or tensor type field ouptut, 
# Type: a string indicates component number or invariant SymbolicConstant
# pick one of the below and enclose it with ''
# possible component values for vectors: 1, 2, 3
# possible component values for tensors: 11, 22, 33, 12, 13, 23
# possible invaraint values: MAGNITUDE, MISES, TRESCA, PRESS, INV3, 
# MAX_PRINCIPAL, MID_PRINCIPAL, MIN_PRINCIPAL, MAX_INPLANE_PRINCIPAL, 
# MIN_INPLANE_PRINCIPAL, OUTOFPLANE_PRINCIPAL
fieldVarComponent = '11'

# Define the instance of model whose result is of interest
# Type: string
instanceName = 'WEB-1'

# Define a potential search range for nodes of the above instance. The poi 
# location should be enclosed in this potential range all the time. It's
# recommonded to define to reduce processing time if model is large. Leave 
# it blank if all nodes of the instance are of interest. Format follows 
# Abaqus nodal path. 
# Type: string 
#potentialNodeRange = '30, 7028:8226:1, 29, 7014:7027:1, 28, 6615:7013:1, 27, 6601:6614:1, 25'        
potentialNodeRange=''
# End of input

import numpy as np
from random import randint
from odbAccess import *
import visualization
from abaqusConstants import *

odb = openOdb(path=odbName, readOnly=True)
instance = odb.rootAssembly.instances[instanceName]
stepRepo = odb.steps

def nLabelStr2IntTuple(nLabelStr):
    """converts an legal nodeList string to a to a tuple consists of integer 
    labels of those nodes
    
        nLabelStr: a string of Abaqus nodeList 
        format of nodeList is same as the one used to create nodal path in 
        Abaqus visualization
        return: a tuple of integer
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
    """create an Abaqus nodeset object, if the nodeset is not created 
    previously
    
        nodeLabels: a tuple consists of labels of nodes
                    or a string
        nodeSetName: a string of nodeSet name, it will be used 
                    as the key of repository
        instance: an Abaqus instance object
        return: an Abaqus nodeSet object
    """
    if nodeSetName not in instance.nodeSets.keys():
        # empty nodeLabesl creates nodeSet of all nodes of a instance
        if nodeLabels == '' or nodeLabels == ' ':
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

def createElementSet(nodeLabels, elementSetName, instance):
    """create an Abaqus odbSet object, if the element set is not created previously
        
        nodeLabels: a tuple consists of labels of elements
        elementSetName: a string of nodeSet name, it will be used 
                        as the key of repository 
        instance: an Abaqus instance object
        return: an Abaqus elementSet object
    """
    if elementSetName not in instance.elementSets.keys():
        elements = ()
        for nLabel in nodeLabels:
            elements += sharedByElementsDict[nLabel] 
        elementSet = instance.ElementSetFromElementLabels(name=elementSetName, elementLabels=elements)
    else:
        elementSet = instance.elementSets[elementSetName]
    return elementSet 

def getNLabelsCoordSeq(frame, nodeSet):
    """ get a tuple of node labels and their corresponded coordiantes
        
        frame: an Abaqus frame object
        nodeSet: an Abaqus nodeSet object
        return: nLabelsSeq: a tuple of node labels
                coordSeq: a tuple of node's coordinates (a 3 element-tuple)
                          in the sequence of node appearing in nLabelSeq
    """
    nLabelsSeq = ()     # in ascending sequence of node label number
    coordSeq = ()       # in cooresponded sequence of nLabelsSeq
    fieldValues = frame.fieldOutputs['U'].getSubset(region=nodeSet).values
    for v, n in zip(fieldValues, nodeSet.nodes):
        # if problem is 2D, displacements is 2D, add 0 for 3rd component
        if len(v.data) == 2:
            coordSeq += (v.data + n.coordinates[0:2], )
        else:
            coordSeq += (v.data + n.coordinates, )
        nLabelsSeq += (n.label, )
    
    return nLabelsSeq, coordSeq

def createCoordField(nodeSet, instance, stepKey):
    """ create an Abaqus scratchOdb that has field output 'COORD' as
    nodes' coordiantes in every frame
    
        nodeSet: an Abaqus nodeSet object
        instance: an Abaqus instance object
        return: None
                created scratchOdb will be accessed through session
    """
    scratchOdb = session.ScratchOdb(odb=odb)
    step = stepRepo[stepKey]
    try:
        ghostStep = scratchOdb.Step(name='Ghost'+stepKey, 
                                    description='Ghost step for your step',
                                    domain=TIME, timePeriod=step.timePeriod)
    except OdbError:
        print 'Ghost step has been created before. Old ghost step will be used.'
        ghostStep = scratchOdb.steps['Ghost'+stepKey]
    if len(ghostStep.frames) != len(step.frames):
        frameNum = 0
        for frame in step.frames:
            ghostFrame = ghostStep.Frame(
                frameId=frame.frameId, 
                frameValue=frame.frameValue,
                description='Ghost frame')
            ghostField = ghostFrame.FieldOutput(
                name='COORD', 
                description='initial coordinates plus displacements',
                type=VECTOR)
            nLabelsSeq, coordSeq = getNLabelsCoordSeq(frame, nodeSet)
            ghostField.addData(position=NODAL, instance=instance, 
                               labels=nLabelsSeq , data=coordSeq)
            frameNum += 1
            print ''.join(('Prepare coordinates for ', stepKey, ': ', 
                           str(frameNum), '/', str(len(step.frames))))
    return None

def closestNodeLabel(nodeSet, poi, step, frame, instance):
    """find the node in nodalRange that has the shortest distance to poi
        
        nodeSet: an Abaqus nodeSet object
        poi: a tuple of coordinates of spatial point of interest
        step: an Abaqus step object
        fraem: an Abaqus frame object
        instance: an Abaqus instance object
        return: an int label of the shortest distance node
    """
    dist = float('inf')
    if 'COORD' in frame.fieldOutputs.keys():
        coordFieldPotential = frame.fieldOutputs['COORD'].getSubset(region=nodeSet)
    else:
        ghostOdbName = [k for k in session.scratchOdbs.keys() if odbName in k][0]
        coordFieldPotential = session.scratchOdbs[ghostOdbName].steps[
            'Ghost'+step.name].getFrame(
            frameValue=frame.frameValue).fieldOutputs['COORD']
    if len(coordFieldPotential.values[0].data) == 2:
        poiA = np.array(poi[0:2])
    else:
        poiA = np.array(poi)
    for nCoordValue in coordFieldPotential.values:
        nCoordData = nCoordValue.data
        distNew = np.linalg.norm(poiA-nCoordData)
        if distNew < dist:
            dist = distNew
            closestNLabel = nCoordValue.nodeLabel
    return closestNLabel

def createXYDataObj(xySequence, xyDataName):
    """ create an Abaqus XYData object
    
        xySequence: a tuple of two tuples of x-data and y-data
        xyDataName: a string indicates XYData name
        return: an Abaqus XYData object
    """
    if xyDataName not in session.xyDataObjects.keys():
        resXYData = session.XYData(data=xySequence, name=xyDataName)
    else:
        session.xyDataObjects[xyDataName].setValues(data=xySequence)
        resXYData = session.xyDataObjects[xyDataName]
    return resXYData

def intString(s):
    """ check if the string s represents a number for components of field
        output
        
        s: a string
        return: True if s represents a component number
                False if s represents invariants
    """
    try: 
        int(s)
        return True
    except ValueError:
        return False
    return None

def sharedByElements(nodeSetPotential, instance):
    """ create a dictionary: node label as keys, a tuple of labels of elements 
    that has the node as values

        nodeSetPotential: an Abaqus nodeSet object
        instance: an Abaqus instance object
        reutrn: a dictionary: key: node labels
                              values: a tuple of labels of elements 
                                      share the node
    """
    nodalElements = {}
    for element in instance.elements:
        for nLabel in element.connectivity:
            if nLabel not in nodalElements.keys():
                nodalElements[nLabel] = ()
            nodalElements[nLabel]+=(element.label,)
    return nodalElements

def fieldVarComponentConstant(fieldVarComponent):
    """ convert the input string-type field component to Abaqus symbolicConstant
        
        fieldComponents: the input string indicates requested 
                         field output component 
        return: an Abaqus symbolicConstant
    """
    if fieldVarComponent == 'MAGNITUDE': return MAGNITUDE
    if fieldVarComponent == 'MISES': return MISES
    if fieldVarComponent == 'TRESCA': return TRESCA
    if fieldVarComponent == 'PRESS': return PRESS
    if fieldVarComponent == 'INV3': return INV3
    if fieldVarComponent ==  'MAX_PRINCIPAL': return MAX_PRINCIPAL
    if fieldVarComponent == 'MAX_PRINCIPAL': return MAX_PRINCIPAL
    if fieldVarComponent == 'MID_PRINCIPAL': return MID_PRINCIPAL
    if fieldVarComponent == 'MIN_PRINCIPAL': return MIN_PRINCIPAL
    if fieldVarComponent == 'MAX_INPLANE_PRINCIPAL': return MAX_INPLANE_PRINCIPAL
    if fieldVarComponent == 'MIN_INPLANE_PRINCIPAL': return MIN_INPLANE_PRINCIPAL
    if fieldVarComponent == 'OUTOFPLANE_PRINCIPAL': return OUTOFPLANE_PRINCIPAL
    return None

def getVarValue(fieldVarName, fieldVarComponent, frame, regionSet):
    """ get the field output value for region in a frame
        
        fieldVarName: the input field output name
        fieldVarComponent: the input field output component
        frame: an Abaqus frame object
        regionSet: an Abaqus nodeSet object or elementSet
        return: the value at node if type of requested field output is NODAL
                the mean value of field output among regionSet 
                    if type of INTEGRATION_POINT
    """
    position = frame.fieldOutputs[fieldVarName].locations[0].position
    if position == NODAL and intString(fieldVarComponent):
        varValue = frame.fieldOutputs[fieldVarName].getSubset(region=regionSet).values[0].data[int(fieldVarComponent)-1]
    elif position == INTEGRATION_POINT and intString(fieldVarComponent):
        fieldVar = ''.join((fieldVarName, fieldVarComponent))
        fieldComponents = frame.fieldOutputs[fieldVarName].getSubset(region=regionSet).getScalarField(componentLabel=fieldVar)
        yDataList = [v.data for v in fieldComponents.values]
        varValue = sum(yDataList)/float(len(yDataList))
    elif not intString(fieldVarComponent):
        fieldComponents = frame.fieldOutputs[fieldVarName].getSubset(region=regionSet).getScalarField(invariant=fieldVarComponentConstant(fieldVarComponent))
        yDataList = [v.data for v in fieldComponents.values]
        varValue = sum(yDataList)/float(len(yDataList))
    return varValue

def plotData(spatialXYData, fieldVarName, fieldVarComponent):
    """ plot data in Abaqus/Visualization
    
        spatialXYData: XYData containing results at spatial location
        xyDataName: plot name
    """
    xyDataName = ''.join((fieldVarName, fieldVarComponent))
    xyPlotName = xyDataName.replace('.', '')
    if xyPlotName not in session.xyPlots.keys():
        xyPlot = session.XYPlot(xyPlotName)
    else:
        xyPlot = session.xyPlots[xyPlotName]
    chartName = session.xyPlots[xyPlotName].charts.keys()[0]
    chart = xyPlot.charts[chartName]
    curve = session.Curve(spatialXYData)
    chart.setValues(curvesToPlot=(curve,))
    
    # set title
    if intString(fieldVarComponent):
        variableName = ''.join((fieldVarName, fieldVarComponent))
    else:
        variableName = ''.join((fieldVarName, '-', fieldVarComponent))
    poiStr = ''.join(('(', str(poi[0]), ', ', str(poi[1]), ', ', str(poi[2]), ')'))
    xyPlot.title.setValues(text=''.join((variableName, ' @ ', poiStr)))
    
    # set axes
    chart.axes1[0].axisData.setValues(title='Time')
    chart.axes2[0].axisData.setValues(title=variableName)
    session.viewports['Viewport: 1'].setValues(displayedObject=xyPlot)
    return None

xySeq = ()
nodeSetName = 'potentialNodeSetForSearch'
nodeSetPotential = createNodeSet(potentialNodeRange, nodeSetName, instance)
fieldOutputInOdb = False    # flag if output field is not in every step
for stepKey in stepRepo.keys():
    step = stepRepo[stepKey]
    if fieldVarName not in step.frames[-1].fieldOutputs.keys():
        print ''.join(('Requested ', fieldVarName, ' is not available in step ',
                        stepKey, ' . Step time skipped.'))
        continue
    else:
        fieldOutputInOdb = True
    
    # create scratchOdb contains coordinates of nodes in the nodeSetPotential
    # at each frame if coordinates are not in field outputs
    if 'COORD' not in step.frames[-1].fieldOutputs.keys():
        createCoordField(nodeSetPotential, instance, stepKey)
    
    # if field output requested is a tensor, create a dictionary with 
    # key=nodeLabel, item=elementsHaveNode
    variablePosition = step.frames[-1].fieldOutputs[fieldVarName].locations[0].position
    if variablePosition == INTEGRATION_POINT:
        sharedByElementsDict = sharedByElements(nodeSetPotential, instance)
    
    frameNum = 0
    stepFrameNum = len(step.frames)
    for frame in step.frames:
        dataPoints = []
        time = step.totalTime + frame.frameValue
        dataPoints.append(time)
        nLabel = closestNodeLabel(nodeSetPotential, poi, step, frame, instance)
        regionSetName = str(nLabel) + stepKey + str(frame.frameId)
        if variablePosition == NODAL:
            regionSet = createNodeSet((nLabel,), regionSetName, instance)
        elif variablePosition == INTEGRATION_POINT:
            regionSet = createElementSet((nLabel,), regionSetName, instance)
        output = getVarValue(fieldVarName, fieldVarComponent, frame, regionSet)
        dataPoints.append(output)
        xySeq += (tuple(dataPoints),)
        frameNum += 1
        print ''.join((stepKey, ': ', str(frameNum), '/', str(stepFrameNum)))
    xyDataName = ''.join((fieldVarName, fieldVarComponent))
    spatialXYData = createXYDataObj(xySequence=xySeq, xyDataName=xyDataName)
    plotData(spatialXYData, fieldVarName, fieldVarComponent)

if not fieldOutputInOdb:
    print ''.join(('Field output ', fieldVarName, ' is not in ODB.'))
