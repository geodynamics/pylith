#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2004  All Rights Reserved
#
#  Copyright 2004 Rensselaer Polytechnic Institute.
#  All worldwide rights reserved.  A license to use, copy, modify and
#  distribute this software for non-commercial research purposes only
#  is hereby granted, provided that this copyright notice and
#  accompanying disclaimer is not modified or removed from the software.
#
#  DISCLAIMER:  The software is distributed "AS IS" without any express
#  or implied warranty, including but not limited to, any implied
#  warranties of merchantability or fitness for a particular purpose
#  or any warranty of non-infringement of any current or pending patent
#  rights.  The authors of the software make no representations about
#  the suitability of this software for any particular purpose.  The
#  entire risk as to the quality and performance of the software is with
#  the user.  Should the software prove defective, the user assumes the
#  cost of all necessary servicing, repair or correction.  In
#  particular, neither Rensselaer Polytechnic Institute, nor the authors
#  of the software are liable for any indirect, special, consequential,
#  or incidental damages related to the software, to the maximum extent
#  the law permits.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from pyre.components.Component import Component


class Pylith3d_scan(Component):


    def __init__(self):
        Component.__init__(self, "pl3dscan", "scanner")

        from pyre.units.pressure import Pa
        from pyre.units.length import m
        from pyre.units.time import s

	print ""
        print "Hello from pl3dscan.__init__ (begin)!"
        print "Setting default keyword values:"

        # default values for extra input (category 2)
        # these can be overriden using a script or keyword=value file

        self.winklerScaleX = 1.0
        self.winklerScaleY = 1.0
        self.winklerScaleZ = 1.0
        self.stressTolerance = 1.0e-12*Pa
        self.minimumStrainPerturbation = 1.0e-7
        self.initialStrainPerturbation = 1.0e-1
        self.preconditionerType = "diagonalNoUpdate"
        self.maxPcgIterations = 3000
        self.displacementAccuracyMult = 1.0
        self.forceAccuracyMult = 1.0
        self.energyAccuracyMult = 1.0
        self.minDisplacementAccuracy = 1.0e-8
        self.minForceAccuracy = 1.0e-8
        self.minEnergyAccuracy = 1.0e-14
        self.quadratureOrder = "Full"
        self.gravityX = 0.0*m/(s*s)
        self.gravityY = 0.0*m/(s*s)
        self.gravityZ = 0.0*m/(s*s)
        self.prestressAutoCompute = False
        self.prestressAutoComputePoisson = -0.49
        self.prestressScaleXx = 1.0
        self.prestressScaleYy = 1.0
        self.prestressScaleZz = 1.0
        self.prestressScaleXy = 1.0
        self.prestressScaleXz = 1.0
        self.prestressScaleYz = 1.0
        self.winklerSlipScaleX = 1.0
        self.winklerSlipScaleY = 1.0
        self.winklerSlipScaleZ = 1.0
        self.f77StandardInput = 5
        self.f77StandardOutput = 6
        self.f77FileInput = 10
        self.f77AsciiOutput = 11
        self.f77PlotOutput = 12

	print ""
        print "Hello from pl3dscan.__init__ (end)!"
        
        return

# derived or automatically-specified quantities (category 3)

    def _init(self, parent):

        from Materials import Materials
        from KeywordValueParse import KeywordValueParse
        import pyre.units
        import pylith3d
        import string
        import os

        uparser = pyre.units.parser()
        matinfo = Materials()
        keyparse = KeywordValueParse()

	print ""
        print "Hello from pl3dscan._init (begin)!"
        print "Scanning ascii files to determine dimensions:"

        # Initialization of all parameters
	# Memory size variable to keep approximate track of all
	# allocated memory.  This does not include python variables and
	# lists.
	self._memorySize = 0L
	self._intSize = 4L
	self._doubleSize = 8L
        # Parameters that are invariant for this geometry type
        self._geometryType = ""
        self._geometryTypeInt = 0
        self._numberSpaceDimensions = 0
        self._numberDegreesFreedom = 0
        self._stateVariableDimension = 0
        self._materialMatrixDimension = 0
        self._numberSkewDimensions = 0
        self._numberSlipDimensions = 0
        self._numberSlipNeighbors = 0
        self._numberTractionDirections = 0
        self._listIddmat = [0]

        # Invariant parameters related to element type
        self._maxElementNodes = 0
        self._maxGaussPoints = 0
        self._maxElementEquations = 0
        self._numberElementTypes = 0
        self._numberElementTypesBase = 0
        self._numberElementNodesBase = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self._pointerToListArrayNumberElementNodesBase = None
        self._pointerToElementTypeInfo = None

        # Invariant parameters related to material model
        self._maxMaterialModels = 0
        self._maxStateVariables = 0
        self._pointerToMaterialModelInfo = None
        self._pointerToMaterialModelStates = None

        # Parameters derived from values in the inventory or the
        # category 2 parameters above.
        self._quadratureOrderInt = 0
        self._analysisTypeInt = 0
        self._prestressAutoComputeInt = 0
        self._pointerToSh = None
        self._pointerToShj = None
        self._pointerToGauss = None

        # Parameters derived from the number of entries in a file

        self._numberNodes = 0
        self._coordinateUnits = "coordinateUnitsInitial12345678"
        self._coordinateScaleFactor = 0.0

        self._numberBcEntries = 0
        self._displacementUnits = "displacementUnitsInitial123456"
        self._displacementScaleFactor = 0.0
        self._velocityUnits = "velocityUnitsInitial1234567890"
        self._velocityScaleFactor = 0.0
        self._forceUnits = "forceUnitsInitial1234567890123"
        self._forceScaleFactor = 0.0

	self._winklerInfo = [0, 0]
        self._numberWinklerEntries = 0
        self._numberWinklerForces = 0

        self._numberRotationEntries = 0
        self._rotationUnits = "rotationUnitsInitial1234567890"
        self._rotationScaleFactor = 0.0

        self._timeStepInfo = [0, 0]
        self._numberTimeStepGroups = 0
        self._totalNumberTimeSteps = 0
        self._timeUnits = "timeUnitsInitial12345678901234"
        self._timeScaleFactor = 0.0

        self._numberFullOutputs = 0

        self._numberLoadHistories = 0

        self._numberMaterials = 0
        self._propertyListSize = 0
        self._propertyList = [0]
        self._pointerToListArrayPropertyList = None
        self._propertyListIndex = [0]
        self._pointerToListArrayPropertyListIndex = None
        self._materialModel = [0]
        self._pointerToListArrayMaterialModel = None
        self._pointerToMaterialInfo = None

        self._elementInfo = [0, 0]
        self._numberElements = 0
        self._connectivitySize = 0

        self._numberPrestressEntries = 0

        self._numberTractionBc = 0
        self._tractionBcUnits = "tractionBcUnitsInitial12345678"
        self._tractionBcScaleFactor = 0.0

        self._numberSplitNodeEntries = 0

        self._numberSlipperyNodeEntries = 0
        self._numberDifferentialForceEntries = 0
	self._slipperyWinklerInfo = [0, 0]
        self._numberSlipperyWinklerEntries = 0
        self._numberSlipperyWinklerForces = 0


        # First see if there is a keyword = value file, which may be used
        # to override parameters from __init__.
        
        if self.inventory.keywordEqualsValueFile == "None":
            self._keywordEqualsValueFile = self.inventory.fileRoot + ".keyval"
        else:
            self._keywordEqualsValueFile = self.inventory.keywordEqualsValueFile

        if os.path.isfile(self._keywordEqualsValueFile):
            file=open(self._keywordEqualsValueFile, 'r')
            while 1:
                line = file.readline()
                if not line: break
                keyvals = keyparse.parseline(line)
                if keyvals[2]:
		    # print keyvals
		    if type(keyvals[1]) == str:
			# print 'Evaluated as string'
		        exec 'self.' + keyvals[0] + '=' + "keyvals[1]"
		    else:
			# print 'Evaluated as other'
		        exec 'self.' + keyvals[0] + '=' + str(keyvals[1])
            file.close()

        # Define information needed from other functions:
        f77FileInput = self.f77FileInput
        prestressAutoCompute = self.prestressAutoCompute
        quadratureOrder = self.quadratureOrder
        
        analysisType = self.inventory.analysisType


        if self.inventory.asciiOutputFile == "None":
            self._asciiOutputFile = self.inventory.fileRoot + ".ascii"
        else:
            self._asciiOutputFile = self.inventory.asciiOutputFile

        if self.inventory.plotOutputFile == "None":
            self._plotOutputFile = self.inventory.fileRoot + ".plot"
        else:
            self._plotOutputFile = self.inventory.plotOutputFile

        if self.inventory.coordinateInputFile == "None":
            self._coordinateInputFile = self.inventory.fileRoot + ".coord"
        else:
            self._coordinateInputFile = self.inventory.coordinateInputFile

        if self.inventory.bcInputFile == "None":
            self._bcInputFile = self.inventory.fileRoot + ".bc"
        else:
            self._bcInputFile = self.inventory.bcInputFile

        if self.inventory.winklerInputFile == "None":
            self._winklerInputFile = self.inventory.fileRoot + ".wink"
        else:
            self._winklerInputFile = self.inventory.winklerInputFile

        if self.inventory.rotationInputFile == "None":
            self._rotationInputFile = self.inventory.fileRoot + ".skew"
        else:
            self._rotationInputFile = self.inventory.rotationInputFile

        if self.inventory.timeStepInputFile == "None":
            self._timeStepInputFile = self.inventory.fileRoot + ".time"
        else:
            self._timeStepInputFile = self.inventory.timeStepInputFile

        if self.inventory.fullOutputInputFile == "None":
            self._fullOutputInputFile = self.inventory.fileRoot + ".fuldat"
        else:
            self._fullOutputInputFile = self.inventory.fullOutputInputFile

        if self.inventory.stateVariableInputFile == "None":
            self._stateVariableInputFile = self.inventory.fileRoot + ".statevar"
        else:
            self._stateVariableInputFile = self.inventory.stateVariableInputFile

        if self.inventory.loadHistoryInputFile == "None":
            self._loadHistoryInputFile = self.inventory.fileRoot + ".hist"
        else:
            self._loadHistoryInputFile = self.inventory.loadHistoryInputFile

        if self.inventory.materialPropertiesInputFile == "None":
            self._materialPropertiesInputFile = self.inventory.fileRoot + ".prop"
        else:
            self._materialPropertiesInputFile = self.inventory.materialPropertiesInputFile

        if self.inventory.materialHistoryInputFile == "None":
            self._materialHistoryInputFile = self.inventory.fileRoot + ".mhist"
        else:
            self._materialHistoryInputFile = self.inventory.materialHistoryInputFile

        if self.inventory.connectivityInputFile == "None":
            self._connectivityInputFile = self.inventory.fileRoot + ".connect"
        else:
            self._connectivityInputFile = self.inventory.connectivityInputFile

        if self.inventory.prestressInputFile == "None":
            self._prestressInputFile = self.inventory.fileRoot + ".prestr"
        else:
            self._prestressInputFile = self.inventory.prestressInputFile

        if self.inventory.tractionInputFile == "None":
            self._tractionInputFile = self.inventory.fileRoot + ".tract"
        else:
            self._tractionInputFile = self.inventory.tractionInputFile

        if self.inventory.splitNodeInputFile == "None":
            self._splitNodeInputFile = self.inventory.fileRoot + ".split"
        else:
            self._splitNodeInputFile = self.inventory.splitNodeInputFile

        if self.inventory.slipperyNodeInputFile == "None":
            self._slipperyNodeInputFile = self.inventory.fileRoot + ".slip"
        else:
            self._slipperyNodeInputFile = self.inventory.slipperyNodeInputFile

        if self.inventory.differentialForceInputFile == "None":
            self._differentialForceInputFile = self.inventory.fileRoot + ".diff"
        else:
            self._differentialForceInputFile = self.inventory.differentialForceInputFile

        if self.inventory.slipperyWinklerInputFile == "None":
            self._slipperyWinklerInputFile = self.inventory.fileRoot + ".winkx"
        else:
            self._slipperyWinklerInputFile = self.inventory.slipperyWinklerInputFile

        # This is a test version where the geometry type is automatically
        # specified by using Pylith3d.  The geometry type is only used for
        # f77 routines and not in pyre. An integer value is also defined
        # for use in f77 routines.
        # Define some integer values that are derived from string variables.

        # Parameters that are invariant for this geometry type
        self._geometryType = "3D"
        self._geometryTypeInt = 4
        self._numberSpaceDimensions = 3
        self._numberDegreesFreedom = 3
        self._stateVariableDimension = 6
        self._materialMatrixDimension = 21
        self._numberSkewDimensions = 2
        self._numberSlipDimensions = 5
        self._numberSlipNeighbors = 4
        self._numberTractionDirections = 2
        # self._listIddmat = [
        #     1, 2, 3, 4, 5, 6,
        #     2, 7, 8, 9,10,11,
        #     3, 8,12,13,14,15,
        #     4, 9,13,16,17,18,
        #     5,10,14,17,19,20,
        #     6,11,15,18,20,21]
        # Changed this to correspond to BLAS packed symmetric matrix format.
        self._listIddmat = [
             1, 2, 4, 7,11,16,
             2, 3, 5, 8,12,17,
             4, 5, 6, 9,13,18,
             7, 8, 9,10,14,19,
            11,12,13,14,15,20,
            16,17,18,19,20,21]

        # Invariant parameters related to element type
        self._maxElementNodes = 20
        self._maxGaussPoints = 27
        self._maxElementEquations = self._numberDegreesFreedom*self._maxElementNodes
        self._numberElementTypes = 62
        self._numberElementTypesBase = 10
        self._numberElementNodesBase = [8, 7, 6, 5, 4, 20, 18, 15, 13, 10]
        self._pointerToListArrayNumberElementNodesBase = pylith3d.intListToArray(
            self._numberElementNodesBase)
	self._memorySize += self._numberElementTypesBase*self._intSize
        self._pointerToElementTypeInfo = pylith3d.allocateInt(
            4*self._numberElementTypes)
	self._memorySize += 4*self._numberElementTypes*self._intSize
        self._pointerToSh = pylith3d.allocateDouble(
            (self._numberSpaceDimensions+1)*
            self._maxElementNodes*
            self._maxGaussPoints*
            self._numberElementTypes)
	self._memorySize += (self._numberSpaceDimensions+1)* \
            self._maxElementNodes* \
            self._maxGaussPoints* \
            self._numberElementTypes* \
	    self._doubleSize
        self._pointerToShj = pylith3d.allocateDouble(
            (self._numberSpaceDimensions+1)*
            self._maxElementNodes*
            self._maxGaussPoints*
            self._numberElementTypes)
	self._memorySize += (self._numberSpaceDimensions+1)* \
            self._maxElementNodes* \
            self._maxGaussPoints* \
            self._numberElementTypes* \
	    self._doubleSize
        self._pointerToGauss = pylith3d.allocateDouble(
            (self._numberSpaceDimensions+1)*
            self._maxGaussPoints*
            self._numberElementTypes)
	self._memorySize += (self._numberSpaceDimensions+1)* \
            self._maxGaussPoints* \
            self._numberElementTypes* \
	    self._doubleSize

        # Invariant parameters related to material model
        self._maxMaterialModels = 20
        self._maxStateVariables = 4
        self._pointerToMaterialModelInfo = pylith3d.allocateInt(
            5*self._maxMaterialModels)
	self._memorySize += 5*self._maxMaterialModels*self._intSize
        self._pointerToMaterialModelStates = pylith3d.allocateInt(
            self._maxStateVariables*self._maxMaterialModels)
	self._memorySize += self._maxStateVariables* \
	    self._maxMaterialModels*self._intSize
        pylith3d.matmod_def(
            self._pointerToMaterialModelInfo,
            self._pointerToMaterialModelStates)

        # Parameters derived from values in the inventory or the
        # category 2 parameters above.
        if analysisType == "dataCheck":
            self._analysisTypeInt = 0
        elif analysisType == "stiffnessFactor":
            self._analysisTypeInt = 1
        elif analysisType == "elasticSolution":
            self._analysisTypeInt = 2
        elif analysisType == "fullSolution":
            self._analysisTypeInt = 3
        else:
            self._analysisTypeInt = 3

        if prestressAutoCompute:
            self._prestressAutoComputeInt = 1
        else:
            self._prestressAutoComputeInt = 0

        if quadratureOrder == "Full":
            self._quadratureOrderInt = 1
        elif quadratureOrder == "Reduced":
            self._quadratureOrderInt = 2
        elif quadratureOrder == "Selective":
            self._quadratureOrderInt = 3
        else:
            self._quadratureOrderInt = 1

        # print "Just before pylith3d.preshape:"
        pylith3d.preshape(
            self._pointerToSh,
            self._pointerToShj,
            self._pointerToGauss,
            self._pointerToElementTypeInfo,
            self._quadratureOrderInt)
        # print "Just after pylith3d.preshape:"

        # Parameters derived from the number of entries in a file.
        # try:
        # print "Just before pylith3d.scan_coords:"
        self._numberNodes = pylith3d.scan_coords(
            f77FileInput,
            self._coordinateUnits,
            self._coordinateInputFile)
        # print "Just after pylith3d.scan_coords:"

        self._coordinateScaleString = \
                                    uparser.parse(string.strip(self._coordinateUnits))
        self._coordinateScaleFactor = self._coordinateScaleString.value

        self._numberBcEntries = pylith3d.scan_bc(
            f77FileInput,
            self._displacementUnits,
            self._velocityUnits,
            self._forceUnits,
            self._bcInputFile)

        self._displacementScaleString = \
                                      uparser.parse(string.strip(self._displacementUnits))
        self._displacementScaleFactor = self._displacementScaleString.value
        self._velocityScaleString = \
                                  uparser.parse(string.strip(self._velocityUnits))
        self._velocityScaleFactor = self._velocityScaleString.value
        self._forceScaleString = \
                               uparser.parse(string.strip(self._forceUnits))
        self._forceScaleFactor = self._forceScaleString.value

        self._winklerInfo = pylith3d.scan_wink(
            f77FileInput,
            self._winklerInputFile)
        self._numberWinklerEntries = self._winklerInfo[0]
        self._numberWinklerForces = self._winklerInfo[1]

        self._numberRotationEntries = pylith3d.scan_skew(
            f77FileInput,
            self._rotationUnits,
            self._rotationInputFile)

        if self._numberRotationEntries != 0:
            self._rotationScaleString = \
                                      uparser.parse(string.strip(self._rotationUnits))
            self._rotationScaleFactor = self._rotationScaleString.value

        self._timeStepInfo = pylith3d.scan_timdat(
            f77FileInput,
            self._timeUnits,
            self._timeStepInputFile)
        self._numberTimeStepGroups = self._timeStepInfo[0]
        self._totalNumberTimeSteps = self._timeStepInfo[1]

        self._timeScaleString = \
                              uparser.parse(string.strip(self._timeUnits))
        self._timeScaleFactor = self._timeScaleString.value

        self._numberFullOutputs = pylith3d.scan_fuldat(
            self._analysisTypeInt,
            self._totalNumberTimeSteps,
            f77FileInput,
            self._fullOutputInputFile)

        self._numberLoadHistories = pylith3d.scan_hist(
            f77FileInput,
            self._loadHistoryInputFile)

        # print "Just before matinfo.readprop"
        self._numberMaterials = matinfo.readprop(self._materialPropertiesInputFile)
        # print "Just after matinfo.readprop"
        # print "numberMaterials = %i" % self._numberMaterials
        self._propertyList = matinfo.propertyList
        self._propertyListIndex = matinfo.propertyIndex
        self._materialModel = matinfo.materialModel
        self._propertyListSize = len(self._propertyList)
        self._pointerToListArrayPropertyList = pylith3d.doubleListToArray(
            self._propertyList)
        self._memorySize += self._propertyListSize*self._doubleSize
        self._pointerToListArrayPropertyListIndex = pylith3d.intListToArray(
            self._propertyListIndex)
        self._memorySize += self._numberMaterials*self._intSize
        self._pointerToListArrayMaterialModel = pylith3d.intListToArray(
            self._materialModel)
        self._memorySize += self._numberMaterials*self._intSize
        self._pointerToMaterialInfo = pylith3d.allocateInt(
            3*self._numberMaterials)
        self._memorySize += 3*self._numberMaterials*self._intSize

        self._elementInfo = pylith3d.scan_connect(
            self._pointerToListArrayNumberElementNodesBase,
            self._pointerToMaterialInfo,
            self._pointerToMaterialModelInfo,
            self._pointerToListArrayMaterialModel,
            self._pointerToListArrayPropertyListIndex,
            self._numberMaterials,
            f77FileInput,
            self._connectivityInputFile)
        self._numberElements = self._elementInfo[0]
        self._connectivitySize = self._elementInfo[1]
        self._pointerToListArrayMaterialModel = None
        self._pointerToListArrayPropertyListIndex = None
        self._memorySize -= 2*self._numberMaterials*self._intSize

        # self._numberPrestressEntries = pylith3d.scan_prestr(
        #     self._stateVariableDimension,
        #     self._numberPrestressGaussPoints,
        #     self._numberElements,
        #     self._prestressAutoComputeInt,
        #     f77FileInput,
        #     self._prestressInputFile)

        # self._numberTractionBc = pylith3d.scan_traction(
        #     self._numberElementNodes,
        #     self._numberTractionDirections,
        #     self._tractionBcUnits,
        #     f77FileInput,
        #     self._tractionInputFile)

        # if self._numberTractionBc != 0:
        #     self._tractionBcScaleString = \
        #                                 1.0*uparser.parse(string.strip(self._tractionBcUnits))
        #     self._tractionBcScaleFactor = \
        #                                 self._tractionBcScaleString/pyre.units.SI.pascal

        self._numberSplitNodeEntries = pylith3d.scan_split(
            f77FileInput,
            self._splitNodeInputFile)

        self._numberSlipperyNodeEntries = pylith3d.scan_slip(
            f77FileInput,
            self._slipperyNodeInputFile)

        self._numberDifferentialForceEntries = pylith3d.scan_diff(
            self._numberSlipperyNodeEntries,
            f77FileInput,
            self._differentialForceInputFile)

        self._slipperyWinklerInfo = pylith3d.scan_winkx(
            self._numberSlipperyNodeEntries,
            f77FileInput,
            self._slipperyWinklerInputFile)
        self._numberSlipperyWinklerEntries = self._slipperyWinklerInfo[0]
        self._numberSlipperyWinklerForces = self._slipperyWinklerInfo[1]

        # except IOError, error:
            # print "Situation:", error
        # except ValueError, error:
            # print "Situation:", error
        # except:
            # print "Exception from Pylith3d_scan!"
                
	print ""
        print "Hello from pl3dscan._init (end)!"

        return


    class Inventory(Component.Inventory):

        import pyre.properties
        import os
#  I don't think this is being used in this section any more.
#        from pyre.units.pressure import pascal


        inventory = [
            pyre.properties.str(
                "title",
                default="Patchtest 1 for linear hex elements, full quadrature."),

            pyre.properties.str(
                "fileRoot",
                default=os.path.join("..","examples","linhex","patchtest","pt1")),

            pyre.properties.str(
                "keywordEqualsValueFile",
                default="None"),

            pyre.properties.str(
                "asciiOutputFile",
                default="None"),

            pyre.properties.str(
                "plotOutputFile",
                default="None"),

            pyre.properties.str(
                "coordinateInputFile",
                default="None"),

            pyre.properties.str(
                "bcInputFile",
                default="None"),

            pyre.properties.str(
                "winklerInputFile",
                default="None"),

            pyre.properties.str(
                "rotationInputFile",
                default="None"),

            pyre.properties.str(
                "timeStepInputFile",
                default="None"),

            pyre.properties.str(
                "fullOutputInputFile",
                default="None"),

            pyre.properties.str(
                "stateVariableInputFile",
                default="None"),

            pyre.properties.str(
                "loadHistoryInputFile",
                default="None"),

            pyre.properties.str(
                "materialPropertiesInputFile",
                default="None"),

            pyre.properties.str(
                "materialHistoryInputFile",
                default="None"),

            pyre.properties.str(
                "connectivityInputFile",
                default="None"),

            pyre.properties.str(
                "prestressInputFile",
                default="None"),

            pyre.properties.str(
                "tractionInputFile",
                default="None"),

            pyre.properties.str(
                "splitNodeInputFile",
                default="None"),

            pyre.properties.str(
                "slipperyNodeInputFile",
                default="None"),

            pyre.properties.str(
                "differentialForceInputFile",
                default="None"),

            pyre.properties.str(
                "slipperyWinklerInputFile",
                default="None"),

            pyre.properties.str(
                "asciiOutput",
                default="full",
                validator=pyre.properties.choice(
                [
                 "none",
                 "echo",
                 "full"])),

            pyre.properties.str(
                "plotOutput",
                default="binary",
                validator=pyre.properties.choice(
                [
                 "ascii",
                 "binary"])),

# Eliminating this option for now, as the geometry type is automatically
# specified by using Pylith3d.
#            pyre.properties.str(
#                "geometryType", default="3D",
#                validator=pyre.properties.choice(["axisymmetric",
#                                                  "planeStrain",
#                                                  "planeStress",
#                                                  "outOfPlane",
#                                                  "3D"])),

            pyre.properties.str(
                "analysisType",
                default="fullSolution",
                validator=pyre.properties.choice(
                [
                 "dataCheck",
                 "stiffnessFactor",
                 "elasticSolution",
                 "fullSolution"])),

            pyre.properties.bool(
                "debuggingOutput",
                default=False),

            pyre.properties.bool(
                "autoRotateSlipperyNodes",
                default=True),

            pyre.properties.int(
                "numberCycles",
                default=1)

            ]


# version
# $Id: Pylith3d_scan.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $

# End of file 
