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
#
# The main function of this code is to emulate the original functionality of
# input.f in the original version of TECTON.  This code segment controls the
# allocation of memory and the reading of the input file.  Additional functions
# covered by this code include the sparse matrix setup portion, which also does
# some memory allocation.  Additional code sections will call the main elastic
# and time-dependent solution drivers, which are presently f77 subroutines.
#
# The code here should be executed after all initializations in Lithomop3d_scan.py
# have been performed, including reading of the keyword=value file.
#


from pyre.components.Component import Component


class Lithomop3d_setup(Component):

    def initialize(self, scanner):
        lm3dscan = scanner

	print ""
        print "Hello from lm3dsetup.initialize (begin)!"
        print "Importing values from scanning phase:"

        # Parameters needed from Lithomop3d_scan.py.
        # These parameters, which may be modified in the run function, should all be
        # available to Lithomop3d_run.py.  Inventory items are not altered here and
        # may be accessed from Lithomop3d_scan.py.

        self.winklerScaleX = lm3dscan.winklerScaleX
        self.winklerScaleY = lm3dscan.winklerScaleY
        self.winklerScaleZ = lm3dscan.winklerScaleZ

        self.stressTolerance = lm3dscan.stressTolerance.value
        self.minimumStrainPerturbation = lm3dscan.minimumStrainPerturbation
        self.initialStrainPerturbation = lm3dscan.initialStrainPerturbation
        self.preconditionerType = lm3dscan.preconditionerType
        self.maxPcgIterations = lm3dscan.maxPcgIterations

        self.displacementAccuracyMult = lm3dscan.displacementAccuracyMult
        self.forceAccuracyMult = lm3dscan.forceAccuracyMult
        self.energyAccuracyMult = lm3dscan.energyAccuracyMult
        
        self.minDisplacementAccuracy = lm3dscan.minDisplacementAccuracy
        self.minForceAccuracy = lm3dscan.minForceAccuracy
        self.minEnergyAccuracy = lm3dscan.minEnergyAccuracy

        self.gravityX = lm3dscan.gravityX.value
        self.gravityY = lm3dscan.gravityY.value
        self.gravityZ = lm3dscan.gravityZ.value

        self.prestressAutoComputePoisson = lm3dscan.prestressAutoComputePoisson
        self.prestressScaleXx = lm3dscan.prestressScaleXx
        self.prestressScaleYy = lm3dscan.prestressScaleYy
        self.prestressScaleZz = lm3dscan.prestressScaleZz
        self.prestressScaleXy = lm3dscan.prestressScaleXy
        self.prestressScaleXz = lm3dscan.prestressScaleXz
        self.prestressScaleYz = lm3dscan.prestressScaleYz

        self.winklerSlipScaleX = lm3dscan.winklerSlipScaleX
        self.winklerSlipScaleY = lm3dscan.winklerSlipScaleY
        self.winklerSlipScaleZ = lm3dscan.winklerSlipScaleZ

        self.f77StandardInput = lm3dscan.f77StandardInput
        self.f77StandardOutput = lm3dscan.f77StandardOutput
        self.f77FileInput = lm3dscan.f77FileInput
        self.f77AsciiOutput = lm3dscan.f77AsciiOutput
        self.f77PlotOutput = lm3dscan.f77PlotOutput
        self.f77UcdOutput = lm3dscan.f77UcdOutput

                                                                   
        # Initialize and define some integer parameters based on string or logical parameters in python

        self.asciiOutputInt = 0
        if lm3dscan.inventory.asciiOutput == "none":
            self.asciiOutputInt = 0
        elif lm3dscan.inventory.asciiOutput == "echo":
            self.asciiOutputInt = 1
        else:
            self.asciiOutputInt = 2
            
        self.plotOutputInt = 0
        if lm3dscan.inventory.plotOutput == "none":
            self.plotOutputInt = 0
        elif lm3dscan.inventory.plotOutput == "ascii":
            self.plotOutputInt = 1
        else:
            self.plotOutputInt = 2
            
        self.ucdOutputInt = 0
        if lm3dscan.inventory.ucdOutput == True:
            self.ucdOutputInt = 1
            
        self.debuggingOutputInt = 0
        if lm3dscan.inventory.debuggingOutput:
            self.debuggingOutputInt = 1
        else:
            self.debuggingOutputInt = 0

        self.autoRotateSlipperyNodesInt = 0
        if lm3dscan.inventory.autoRotateSlipperyNodes:
            self.autoRotateSlipperyNodesInt = 2
        else:
            self.autoRotateSlipperyNodesInt = 1

        self.preconditionerTypeInt = 0
        if lm3dscan.preconditionerType == "diagonalNoUpdate":
            self.preconditionerTypeInt = 1
        elif lm3dscan.preconditionerType == "gaussSeidelNoUpdate":
            self.preconditionerTypeInt = 2
        elif lm3dscan.preconditionerType == "diagonalUpdate":
            self.preconditionerTypeInt = 3
        elif lm3dscan.preconditionerType == "gaussSeidelUpdate":
            self.preconditionerTypeInt = 4
        else:
            self.preconditionerTypeInt = 1

        # Get some parameters from the inventory list.

        self.title = lm3dscan.inventory.title
        self.numberCycles = lm3dscan.inventory.numberCycles

        # Category 2 parameters needed from lm3dscan._init.  These parameters are
        # not meant to be altered by the user.
        
	# Approximate memory allocation info:
	self.memorySize = lm3dscan._memorySize
	self.intSize = lm3dscan._intSize
	self.doubleSize = lm3dscan._doubleSize
        # First get all filenames

        self.asciiOutputFile = lm3dscan._asciiOutputFile
        self.plotOutputFile = lm3dscan._plotOutputFile
        self.ucdOutputRoot = lm3dscan._ucdOutputRoot
        self.coordinateInputFile = lm3dscan._coordinateInputFile
        self.bcInputFile = lm3dscan._bcInputFile
        self.winklerInputFile = lm3dscan._winklerInputFile
        self.rotationInputFile = lm3dscan._rotationInputFile
        self.timeStepInputFile = lm3dscan._timeStepInputFile
        self.fullOutputInputFile = lm3dscan._fullOutputInputFile
        self.stateVariableInputFile = lm3dscan._stateVariableInputFile
        self.loadHistoryInputFile = lm3dscan._loadHistoryInputFile
        self.materialPropertiesInputFile = lm3dscan._materialPropertiesInputFile
        self.materialHistoryInputFile = lm3dscan._materialHistoryInputFile
        self.connectivityInputFile = lm3dscan._connectivityInputFile
        self.prestressInputFile = lm3dscan._prestressInputFile
        self.tractionInputFile = lm3dscan._tractionInputFile
        self.splitNodeInputFile = lm3dscan._splitNodeInputFile
        self.slipperyNodeInputFile = lm3dscan._slipperyNodeInputFile
        self.differentialForceInputFile = lm3dscan._differentialForceInputFile
        self.slipperyWinklerInputFile = lm3dscan._slipperyWinklerInputFile

        # Parameters that are invariant for this geometry type
        self.geometryType = lm3dscan._geometryType
        self.geometryTypeInt = lm3dscan._geometryTypeInt
        self.numberSpaceDimensions = lm3dscan._numberSpaceDimensions
        self.numberDegreesFreedom = lm3dscan._numberDegreesFreedom
        self.stateVariableDimension = lm3dscan._stateVariableDimension
        self.materialMatrixDimension = lm3dscan._materialMatrixDimension
        self.numberSkewDimensions = lm3dscan._numberSkewDimensions
        self.numberSlipDimensions = lm3dscan._numberSlipDimensions
        self.numberSlipNeighbors = lm3dscan._numberSlipNeighbors
        self.numberTractionDirections = lm3dscan._numberTractionDirections
        self.listIddmat = lm3dscan._listIddmat

        # Invariant parameters related to element type
        self.maxElementNodes = lm3dscan._maxElementNodes
        self.maxGaussPoints = lm3dscan._maxGaussPoints
        self.maxElementEquations = lm3dscan._maxElementEquations
        self.numberElementTypes = lm3dscan._numberElementTypes
        self.pointerToListArrayNumberElementNodesBase = lm3dscan._pointerToListArrayNumberElementNodesBase
        self.pointerToElementTypeInfo = lm3dscan._pointerToElementTypeInfo
        self.pointerToSh = lm3dscan._pointerToSh
        self.pointerToShj = lm3dscan._pointerToShj
        self.pointerToGauss = lm3dscan._pointerToGauss

        # Invariant parameters related to material model
        self.maxMaterialModels = lm3dscan._maxMaterialModels
        self.maxStateVariables = lm3dscan._maxStateVariables
        self.pointerToMaterialModelInfo = lm3dscan._pointerToMaterialModelInfo
        self.pointerToMaterialModelStates = lm3dscan._pointerToMaterialModelStates

        # Parameters derived from values in the inventory or the
        # category 2 parameters.
        self.analysisTypeInt = lm3dscan._analysisTypeInt
        self.quadratureOrderInt = lm3dscan._quadratureOrderInt
        self.prestressAutoComputeInt = lm3dscan._prestressAutoComputeInt


        # Parameters derived from the number of entries in a file

        self.numberNodes = lm3dscan._numberNodes
        self.coordinateScaleFactor = lm3dscan._coordinateScaleFactor

        self.numberBcEntries = lm3dscan._numberBcEntries
        self.displacementScaleFactor = lm3dscan._displacementScaleFactor
        self.velocityScaleFactor = lm3dscan._velocityScaleFactor
        self.forceScaleFactor = lm3dscan._forceScaleFactor

        self.numberWinklerEntries = lm3dscan._numberWinklerEntries
        self.numberWinklerForces = lm3dscan._numberWinklerForces

        self.numberRotationEntries = lm3dscan._numberRotationEntries
        self.rotationScaleFactor = lm3dscan._rotationScaleFactor

        self.numberTimeStepGroups = lm3dscan._numberTimeStepGroups
        self.totalNumberTimeSteps = lm3dscan._totalNumberTimeSteps
        self.timeScaleFactor = lm3dscan._timeScaleFactor

        self.numberFullOutputs = lm3dscan._numberFullOutputs

        self.numberLoadHistories = lm3dscan._numberLoadHistories

        self.numberMaterials = lm3dscan._numberMaterials
        self.propertyListSize = lm3dscan._propertyListSize
        self.pointerToListArrayPropertyList = lm3dscan._pointerToListArrayPropertyList
        self.pointerToMaterialInfo = lm3dscan._pointerToMaterialInfo

        self.numberElements = lm3dscan._numberElements
        self.connectivitySize = lm3dscan._connectivitySize

        self.numberPrestressEntries = lm3dscan._numberPrestressEntries

        self.numberTractionBc = lm3dscan._numberTractionBc
        self.tractionBcScaleFactor = lm3dscan._tractionBcScaleFactor

        self.numberSplitNodeEntries = lm3dscan._numberSplitNodeEntries

        self.numberSlipperyNodeEntries = lm3dscan._numberSlipperyNodeEntries
        self.numberDifferentialForceEntries = lm3dscan._numberDifferentialForceEntries
        self.numberSlipperyWinklerEntries = lm3dscan._numberSlipperyWinklerEntries
        self.numberSlipperyWinklerForces = lm3dscan._numberSlipperyWinklerForces

	print ""
        print "Hello from lm3dsetup.initialize (end)!"

        return

    def run(self):
        import lithomop3d

        # This function performs a lot of memory allocation, and also bundles several
        # scalar values into lists.  It is assumed that the scalar values contained in
        # the lists will no longer be accessible (or useful) from within python once
        # they are passed as arrays.  The contents of the lists should be maintained as
        # global variables, however.

	print ""
        print "Hello from lm3dsetup.run (begin)!"
        print "Reading problem definition and allocating storage:"

        # Initialize parameters that are defined in this function

        self.numberGlobalEquations = 0
        self.currentNumberEquations = 0

        self.totalNumberSplitNodes = 0
        self.totalNumberSlipperyNodes = 0

        self.pointerToAlnz = None
        self.pointerToJa = None
        self.pointerToPcg = None
        self.pointerToZcg = None

        self.pointerToB = None
        self.pointerToBtot = None
        self.pointerToBres = None
        self.pointerToPvec = None
        self.pointerToGvec1 = None
        self.pointerToGvec2 = None
        self.listGrav = [0.0, 0.0, 0.0]
        self.pointerToListArrayGrav = None

        self.pointerToX = None
        self.pointerToD = None
        self.pointerToDeld = None
        self.pointerToDprev = None
        self.pointerToDcur = None
        self.pointerToId = None
        self.pointerToIwink = None
        self.pointerToWink = None
        self.listNsysdat = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.pointerToListArrayNsysdat = None

        self.pointerToIbond = None
        self.pointerToBond = None

        self.pointerToDx = None
        self.pointerToDeldx = None
        self.pointerToDxcur = None
        self.pointerToDiforc = None
        self.pointerToIdx = None
        self.pointerToIwinkx = None
        self.pointerToWinkx = None
        self.pointerToIdslp = None
        self.pointerToIpslp = None
        self.pointerToIdhist = None

        self.pointerToFault = None
        self.pointerToNfault = None
        self.pointerToDfault = None
        self.pointerToTfault = None
        self.pointerToIdftn = None

        self.pointerToS = None
        self.pointerToStemp = None

        self.pointerToState = None
        self.pointerToDstate = None
        self.stateSize = 0
        self.pointerToDmat = None
        self.dmatSize = 0
        self.pointerToIen = None
        self.pointerToLm = None
        self.pointerToLmx = None
        self.pointerToLmf = None
        self.pointerToInfiel = None
        self.pointerToListArrayIddmat = None
        self.listNpar = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.pointerToListArrayNpar = None
        self.elementSizeInfo = [0, 0]

	#  The following pointers are temporarily being set to numeric
	#  values to avoid error messages from python.
	#  They are not being used at present, and whend they retain
	#  the python value 'None', they cause problems.
	#  Another option:  Set these to a value using an allocation of
	#  zero length.
        self.pointerToIelno = None
        self.pointerToIside = None
        self.pointerToIhistry = None
        self.pointerToPres = None
        self.pointerToPdir = None

        self.pointerToMhist = None

        self.pointerToHistry = None
        self.listRtimdat = [0.0, 0.0, 0.0]
        self.pointerToListArrayRtimdat = None
        self.listNvisdat = [0, 0, 0, 0]
        self.pointerToListArrayNvisdat = None
        self.pointerToMaxstp = None
        self.pointerToDelt = None
        self.pointerToAlfa = None
        self.pointerToMaxit = None
        self.pointerToNtdinit = None
        self.pointerToLgdef = None
        self.pointerToUtol = None
        self.pointerToFtol = None
        self.pointerToEtol = None
        self.pointerToItmax = None

        self.listRgiter = [0.0, 0.0, 0.0]
        self.pointerToListArrayRgiter = None
        self.listRmin = [0.0, 0.0, 0.0]
        self.pointerToListArrayRmin = None
        self.listRmult = [0.0, 0.0, 0.0]
        self.pointerToListArrayRmult = None
        self.listNsiter = [0, 0]
        self.pointerToListArrayNsiter = None
        
        self.pointerToSkew = None

        self.pointerToIprint = None
        self.listNcodat = [0, 0]
        self.pointerToListArrayNcodat = None
        self.listNunits = [0, 0, 0, 0, 0]
        self.pointerToListArrayNunits = None
        self.listNprint = [0, 0, 0]
        self.pointerToListArrayNprint = None
        self.pointerToIstatout = None


        # Arrays that can be deallocated after use.
        self.pointerToTimes = None
        self.pointerToIndmat = None
        self.pointerToImgrp = None
        self.pointerToNslip = None
        self.pointerToXtmp = None
        self.pointerToItmp = None
        self.pointerToItmp1 = None
        self.pointerToItmp2 = None
        self.pointerToIndx = None
        self.pointerToLink = None
        self.pointerToNbrs = None
        self.listWscal = [0.0, 0.0, 0.0]
        self.pointerToListArrayWscal = None
        self.listPrscal = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.pointerToListArrayPrscal = None
        self.listWxscal = [0.0, 0.0, 0.0]
        self.pointerToListArrayWxscal = None

        # Sparse matrix info
        self.workingArraySize = 0
        self.stiffnessMatrixSize = 0
        self.stiffnessTrueSize = 0
        self.stiffnessOffDiagonalSize = 0
        self.stiffnessMatrixInfo = [0, 0]
        self.minimumNonzeroTermsPerRow = 0
        self.maximumNonzeroTermsPerRow = 0
        self.averageNonzeroTermsPerRow = 0.0
        self.stiffnessMatrixStats = [0, 0, 0.0]


        # Make lists that are used as arrays in the f77 function calls below.

        self.listWscal = [
            self.winklerScaleX,
            self.winklerScaleY,
            self.winklerScaleZ]
        self.pointerToListArrayWscal = lithomop3d.doubleListToArray(
            self.listWscal)
	self.memorySize += 3*self.doubleSize

        self.listRmult = [
            self.displacementAccuracyMult,
            self.forceAccuracyMult,
            self.energyAccuracyMult]
        self.pointerToListArrayRmult = lithomop3d.doubleListToArray(
            self.listRmult)
	self.memorySize += 3*self.doubleSize

        self.listRmin = [
            self.minDisplacementAccuracy,
            self.minForceAccuracy,
            self.minEnergyAccuracy]
        self.pointerToListArrayRmin = lithomop3d.doubleListToArray(
            self.listRmin)
	self.memorySize += 3*self.doubleSize

        self.listGrav = [
            self.gravityX,
            self.gravityY,
            self.gravityZ]
        self.pointerToListArrayGrav = lithomop3d.doubleListToArray(
            self.listGrav)
	self.memorySize += 3*self.doubleSize

        self.listPrscal = [
            self.prestressScaleXx,
            self.prestressScaleYy,
            self.prestressScaleZz,
            self.prestressScaleXy,
            self.prestressScaleXz,
            self.prestressScaleYz]
        self.pointerToListArrayPrscal = lithomop3d.doubleListToArray(
            self.listPrscal)
	self.memorySize += 6*self.doubleSize
                                  
        self.listWxscal = [
            self.winklerSlipScaleX,
            self.winklerSlipScaleY,
            self.winklerSlipScaleZ]
        self.pointerToListArrayWxscal = lithomop3d.doubleListToArray(
            self.listWxscal)
	self.memorySize += 3*self.doubleSize

        # Write out global parameters
        
        lithomop3d.write_global_info(
            self.title,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.numberNodes,
            self.analysisTypeInt,
            self.debuggingOutputInt,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputFile,
            self.plotOutputFile)

        # Node-based info (coordinates, ID arrays, displacement arrays, winkler arrays,
        #  BC, and skew BC).

        self.pointerToX = lithomop3d.allocateDouble(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes* \
	    self.doubleSize
        self.pointerToId = lithomop3d.allocateInt(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes* \
	    self.intSize
        self.pointerToIwink = lithomop3d.allocateInt(
            2*self.numberWinklerForces)
	self.memorySize += 2*self.numberWinklerForces*self.intSize
        self.pointerToWink = lithomop3d.allocateDouble(
            self.numberWinklerForces)
	self.memorySize += self.numberWinklerForces*self.doubleSize
        self.pointerToIbond = lithomop3d.allocateInt(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes* \
	    self.intSize
        self.pointerToBond = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes* \
	    self.doubleSize
        # For now, try allocating everything whether it is used or not.
        #        if self.numberRotationEntries != 0 or self.autoRotateSlipperyNodes:
        #            self.pointerToSkew = lithomop3d.allocateDouble(
        #                self.numberSkewDimensions*self.numberNodes)
        self.pointerToSkew = lithomop3d.allocateDouble(
            self.numberSkewDimensions*self.numberNodes)
	self.memorySize += self.numberSkewDimensions* \
	    self.numberNodes* \
	    self.doubleSize


        # try:

        lithomop3d.read_coords(
            self.pointerToX,
            self.coordinateScaleFactor,
            self.numberNodes,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.coordinateInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)

        self.numberGlobalEquations = lithomop3d.read_bc(
            self.pointerToBond,
            self.displacementScaleFactor,
            self.velocityScaleFactor,
            self.forceScaleFactor,
            self.pointerToIbond,
            self.pointerToId,
            self.numberNodes,
            self.numberBcEntries,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.bcInputFile,
            self.asciiOutputFile)
        self.currentNumberEquations = self.numberGlobalEquations

        lithomop3d.read_wink(
            self.pointerToWink,
            self.pointerToListArrayWscal,
            self.pointerToIwink,
            self.pointerToId,
            self.numberNodes,
            self.numberWinklerForces,
            self.numberWinklerEntries,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.winklerInputFile,
            self.asciiOutputFile)

        lithomop3d.read_skew(
            self.pointerToSkew,
            self.rotationScaleFactor,
            self.numberRotationEntries,
            self.numberNodes,
            self.autoRotateSlipperyNodesInt,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.rotationInputFile,
            self.asciiOutputFile)

        # except IOError, error:
            # print "Situation:", error
        # except ValueError, error:
            # print "Situation:", error
        # except ArithmeticError, error:
            # print "Situation:", error
        # except MemoryError, error:
            # print "Situation:", error
        # except Exception, error:
            # print "Exception in block between read_coords and read_skew!", error

        # Output stress computation and subiteration parameters.
        lithomop3d.write_strscomp(
            self.stressTolerance,
            self.minimumStrainPerturbation,
            self.initialStrainPerturbation,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)
        
        lithomop3d.write_subiter(
            self.pointerToListArrayRmult,
            self.pointerToListArrayRmin,
            self.preconditionerTypeInt,
            self.maxPcgIterations,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)
                             
        # Allocate and read time step, time output, and load history info.
        self.pointerToHistry = lithomop3d.allocateDouble(
            (self.totalNumberTimeSteps+1)*self.numberLoadHistories)
	self.memorySize += (self.totalNumberTimeSteps+1)* \
	    self.numberLoadHistories* \
	    self.doubleSize
        self.pointerToMaxstp = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.intSize
        self.pointerToDelt = lithomop3d.allocateDouble(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.doubleSize
        self.pointerToAlfa = lithomop3d.allocateDouble(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.doubleSize
        self.pointerToMaxit = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.intSize
        self.pointerToNtdinit = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.intSize
        self.pointerToLgdef = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.intSize
        self.pointerToUtol = lithomop3d.allocateDouble(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.doubleSize
        self.pointerToFtol = lithomop3d.allocateDouble(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.doubleSize
        self.pointerToEtol = lithomop3d.allocateDouble(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.doubleSize
        self.pointerToItmax = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.intSize
        self.pointerToIprint = lithomop3d.allocateInt(
            self.numberFullOutputs)
	self.memorySize += self.numberFullOutputs*self.intSize
        self.pointerToTimes = lithomop3d.allocateDouble(
            self.totalNumberTimeSteps+1)
	self.memorySize += (self.totalNumberTimeSteps+1)*self.doubleSize
        self.pointerToIstatout = lithomop3d.allocateInt(
            2*self.maxStateVariables)
	self.memorySize += 4*self.maxStateVariables*self.intSize

        # try:
        lithomop3d.read_timdat(
            self.pointerToDelt,
            self.pointerToAlfa,
            self.pointerToUtol,
            self.pointerToFtol,
            self.pointerToEtol,
            self.pointerToTimes,
            self.timeScaleFactor,
            self.pointerToMaxstp,
            self.pointerToMaxit,
            self.pointerToNtdinit,
            self.pointerToLgdef,
            self.pointerToItmax,
            self.numberTimeStepGroups,
            self.totalNumberTimeSteps,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.timeStepInputFile,
            self.asciiOutputFile)

        lithomop3d.read_fuldat(
            self.pointerToIprint,
            self.numberFullOutputs,
            self.analysisTypeInt,
            self.numberCycles,
            self.totalNumberTimeSteps,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.fullOutputInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)

        lithomop3d.read_stateout(
            self.pointerToIstatout,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.stateVariableInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)

        lithomop3d.read_hist(
            self.pointerToHistry,
            self.pointerToTimes,
            self.numberLoadHistories,
            self.totalNumberTimeSteps,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.loadHistoryInputFile,
            self.asciiOutputFile)
        self.pointerToTimes = None
        self.memorySize -= (self.totalNumberTimeSteps+1)*self.doubleSize
            
        # except IOError, error:
            # print "Situation:", error
        # except ValueError, error:
            # print "Situation:", error
        # except ArithmeticError, error:
            # print "Situation:", error
        # except MemoryError, error:
            # print "Situation:", error
        # except Exception, error:
            # print "Exception in block between read_timdat and read_hist!", error

        # Allocate and read info on material properties, connectivities, and prestresses
        lithomop3d.write_element_info(
            self.numberElements,
            self.quadratureOrderInt,
            self.prestressAutoComputeInt,
            self.prestressAutoComputePoisson,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        self.pointerToIen = lithomop3d.allocateInt(
            self.connectivitySize)
	self.memorySize += self.connectivitySize*self.intSize
        self.pointerToInfiel = lithomop3d.allocateInt(
            6*self.numberElements)
	self.memorySize += 6*self.numberElements*self.intSize
        self.pointerToIndmat = lithomop3d.allocateInt(
            self.numberMaterials)
	self.memorySize += self.numberMaterials*self.intSize
        self.pointerToImgrp = lithomop3d.allocateInt(
            self.numberMaterials)
	self.memorySize += self.numberMaterials*self.intSize
        self.pointerToMhist = lithomop3d.allocateInt(
            self.propertyListSize)
	self.memorySize += self.propertyListSize*self.intSize

        # try:
        lithomop3d.write_props(
            self.pointerToListArrayPropertyList,
            self.pointerToListArrayGrav,
            self.pointerToMaterialInfo,
            self.pointerToMaterialModelInfo,
            self.numberMaterials,
            self.propertyListSize,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputFile,
            self.plotOutputFile)

        # print "Just before read_mathist:"
        # print "memorySize: %d" % self.memorySize
        lithomop3d.read_mathist(
            self.pointerToMhist,
            self.pointerToMaterialInfo,
            self.numberMaterials,
            self.propertyListSize,
            self.numberLoadHistories,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.materialHistoryInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)
        # print "Just before read_connect:"
        # print "pointerToListArrayNumberElementNodesBase: %i" % self.pointerToListArrayNumberElementNodesBase
        # print "pointerToElementTypeInfo: %i" % self.pointerToElementTypeInfo
        # print "pointerToMaterialInfo: %i" % self.pointerToMaterialInfo
        # print "pointerToMaterialModelInfo: %i" % self.pointerToMaterialModelInfo
        # print "pointerToIen: %i" % self.pointerToIen
        # print "pointerToInfiel: %i" % self.pointerToInfiel
        # print "pointerToIndmat: %i" % self.pointerToIndmat
        # print "pointerToImgrp: %i" % self.pointerToImgrp
        # print "connectivitySize: %i" % self.connectivitySize
        # print "numberElements: %i" % self.numberElements
        # print "numberNodes: %i" % self.numberNodes
        # print "numberMaterials: %i" % self.numberMaterials
        # print "f77FileInput: %i" % self.f77FileInput


        self.elementSizeInfo = lithomop3d.read_connect(
            self.pointerToListArrayNumberElementNodesBase,
            self.pointerToElementTypeInfo,
            self.pointerToMaterialInfo,
            self.pointerToMaterialModelInfo,
            self.pointerToIen,
            self.pointerToInfiel,
            self.pointerToIndmat,
            self.pointerToImgrp,
            self.connectivitySize,
            self.numberElements,
            self.numberNodes,
            self.numberMaterials,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.connectivityInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)
        self.stateSize = self.elementSizeInfo[0]
        self.dmatSize = self.elementSizeInfo[1]
        # print "stateSize: %i" % self.stateSize
        # print "dmatSize: %i" % self.dmatSize
        self.pointerToIndmat = None
        self.memorySize -= self.numberMaterials*self.intSize
        self.pointerToImgrp = None
        self.memorySize -= self.numberMaterials*self.intSize

        # write mesh info to UCD file, if requested
        if self.ucdOutputInt == 1:
            lithomop3d.write_ucd_mesh(
                self.pointerToX,
                self.numberNodes,
                self.pointerToIen,
                self.pointerToInfiel,
                self.numberElements,
                self.connectivitySize,
                self.pointerToSh,
                self.pointerToElementTypeInfo,
                self.pointerToIstatout,
                self.f77UcdOutput,
                self.ucdOutputRoot)
                
        # print "Just after read_connect"
        # print "memorySize: %d" % self.memorySize
        # lithomop3d.read_prestr(
        #     self.pointerToStn,
        #     self.pointerToSt0,
        #     self.pointerToListArrayPrscal,
        #     self.numberStressComponents,
        #     self.numberGaussPoints,
        #     self.numberPrestressGaussPoints,
        #     self.numberElements,
        #     self.numberPrestressEntries,
        #     self.prestressAutoComputeInt,
        #     self.asciiOutputInt,
        #     self.f77FileInput,
        #     self.f77AsciiOutput,
        #     self.prestressInputFile,
        #     self.asciiOutputFile)
        # print "Just after read_prestr"

        # except IOError, error:
            # print "Situation:", error
        # except ValueError, error:
            # print "Situation:", error
        # except ArithmeticError, error:
            # print "Situation:", error
        # except MemoryError, error:
            # print "Situation:", error
        # except Exception, error:
            # print "Exception in block between read_prop and read_connect!", error

        # Read traction, split node, and slippery node input files.
        self.pointerToIelno = lithomop3d.allocateInt(
            self.numberTractionBc)
        self.pointerToIside = lithomop3d.allocateInt(
            self.numberTractionBc)
        self.pointerToIhistry = lithomop3d.allocateInt(
            self.numberTractionBc)
	#  Note that the following dimension is definitely wrong.
        self.pointerToPres = lithomop3d.allocateDouble(
            self.numberTractionBc)
        self.pointerToPdir = lithomop3d.allocateDouble(
            self.numberTractionDirections*self.numberTractionBc)
        self.pointerToNfault = lithomop3d.allocateInt(
            3*self.numberSplitNodeEntries)
	self.memorySize += 3*self.numberSplitNodeEntries*self.intSize
        self.pointerToFault = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSplitNodeEntries)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberSplitNodeEntries*self.doubleSize
        self.pointerToNslip = lithomop3d.allocateInt(
            self.numberSlipDimensions*self.numberSlipperyNodeEntries)
	self.memorySize += self.numberSlipDimensions* \
	    self.numberSlipperyNodeEntries*self.intSize
        self.pointerToIdhist = lithomop3d.allocateInt(
            self.numberNodes)
	self.memorySize += self.numberNodes*self.intSize
        self.pointerToDiforc = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        # For now, try allocating everything whether it is used or not.
        #        if self.numberDifferentialForceEntries != 0:
        #            self.pointerToIdhist = lithomop3d.allocateInt(
        #                self.numberNodes)
        #            self.pointerToDiforc = lithomop3d.allocateDouble(
        #                self.numberDegreesFreedom*self.numberNodes)
        # try:
        # lithomop3d.read_traction(
        #     self.pointerToPres,
        #     self.pointerToPdir,
        #     self.tractionBcScaleFactor,
        #     self.pointerToIelno,
        #     self.pointerToIside,
        #     self.pointerToIhistry,
        #     self.numberTractionBc,
        #     self.numberElementNodes,
        #     self.numberTractionDirections,
        #     self.f77FileInput,
        #     self.f77AsciiOutput,
        #     self.asciiOutputInt,
        #     self.tractionInputFile,
        #     self.asciiOutputFile)
        # print "Just after read_traction"

        self.totalNumberSplitNodes = lithomop3d.read_split(
            self.pointerToFault,
            self.pointerToNfault,
            self.numberSplitNodeEntries,
            self.numberNodes,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.splitNodeInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)
        # print "Just after read_split"
        # print "memorySize: %d" % self.memorySize

        self.totalNumberSlipperyNodes = lithomop3d.read_slip(
            self.pointerToNslip,
            self.numberSlipperyNodeEntries,
            self.numberNodes,
            self.autoRotateSlipperyNodesInt,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.slipperyNodeInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)
        # print "Just after read_slip"

        lithomop3d.read_diff(
            self.pointerToDiforc,
            self.pointerToNslip,
            self.pointerToIdhist,
            self.numberSlipperyNodeEntries,
            self.numberDifferentialForceEntries,
            self.numberNodes,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.differentialForceInputFile,
            self.asciiOutputFile)
        # print "Just after read_diff"

        # except IOError, error:
            # print "Situation:", error
        # except ValueError, error:
            # print "Situation:", error
        # except ArithmeticError, error:
            # print "Situation:", error
        # except MemoryError, error:
            # print "Situation:", error
        # except Exception, error:
            # print "Exception in block between read_traction and read_diff!", error

        # Allocate id arrays for split and slippery nodes and differential winkler arrays, then
        # adjust equations for the presence of slippery nodes, and read in the differential
        # winkler forces.
        # Note:  all input is done after this section.  The remainder of this function allocates
        # memory for computations.
        
        self.pointerToIdx = lithomop3d.allocateInt(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes*self.intSize
        self.pointerToIdftn = lithomop3d.allocateInt(
            self.totalNumberSplitNodes)
	self.memorySize += self.totalNumberSplitNodes*self.intSize
        self.pointerToWinkx = lithomop3d.allocateDouble(
            self.numberSlipperyWinklerForces)
	self.memorySize += self.numberSlipperyWinklerForces*self.doubleSize
        self.pointerToIwinkx = lithomop3d.allocateInt(
            2*self.numberSlipperyWinklerForces)
	self.memorySize += 2*self.numberSlipperyWinklerForces*self.intSize
        self.pointerToIdslp = lithomop3d.allocateInt(
            self.numberNodes)
	self.memorySize += self.numberNodes*self.intSize
        self.pointerToIpslp = lithomop3d.allocateInt(
            self.numberSlipNeighbors*self.totalNumberSlipperyNodes)
	self.memorySize += self.numberSlipNeighbors* \
	    self.totalNumberSlipperyNodes*self.intSize

        lithomop3d.id_split(
            self.pointerToNfault,
            self.pointerToIdftn,
            self.numberNodes,
            self.numberSplitNodeEntries,
            self.totalNumberSplitNodes,
            self.f77PlotOutput,
            self.plotOutputInt,
            self.plotOutputFile)
        # print "Just after id_split"
	# print "memorySize: %d" % self.memorySize
        # Deallocating this array since it does not appear to be used elsewhere.
        # It may be useful in the future, though.
        self.pointerToIdftn = None
	self.memorySize -= self.totalNumberSplitNodes*self.intSize
                        
        self.numberGlobalEquations = lithomop3d.adjid(
            self.pointerToId,
            self.pointerToIdx,
            self.pointerToNslip,
            self.pointerToIdslp,
            self.numberSlipperyNodeEntries,
            self.numberNodes,
            self.totalNumberSlipperyNodes,
            self.currentNumberEquations)
        # print "Just after adjid"

        if self.totalNumberSlipperyNodes != 0 and self.autoRotateSlipperyNodesInt ==  2:
            self.pointerToXtmp = lithomop3d.allocateDouble(
                self.totalNumberSlipperyNodes)
	    self.memorySize += self.totalNumberSlipperyNodes*self.doubleSize
            self.pointerToItmp = lithomop3d.allocateInt(
                self.totalNumberSlipperyNodes)
	    self.memorySize += self.totalNumberSlipperyNodes*self.intSize
            self.pointerToItmp1 = lithomop3d.allocateInt(
                self.totalNumberSlipperyNodes)
	    self.memorySize += self.totalNumberSlipperyNodes*self.intSize
            self.pointerToItmp2 = lithomop3d.allocateInt(
                self.totalNumberSlipperyNodes)
	    self.memorySize += self.totalNumberSlipperyNodes*self.intSize

            lithomop3d.nfind(
                self.pointerToX,
                self.pointerToXtmp,
                self.pointerToIdslp,
                self.pointerToIpslp,
                self.pointerToItmp,
                self.pointerToItmp1,
                self.pointerToItmp2,
                self.pointerToNslip,
                self.numberSlipperyNodeEntries,
                self.totalNumberSlipperyNodes,
                self.numberNodes)
            # print "Just after nfind"
	    # print "memorySize: %d" % self.memorySize

            self.pointerToXtmp = None
            self.pointerToItmp = None
            self.pointerToItmp1 = None
            self.pointerToItmp2 = None
	    self.memorySize -= self.totalNumberSlipperyNodes*self.doubleSize
	    self.memorySize -= self.totalNumberSlipperyNodes*self.intSize
	    self.memorySize -= self.totalNumberSlipperyNodes*self.intSize
	    self.memorySize -= self.totalNumberSlipperyNodes*self.intSize
                         
        # try:
        lithomop3d.read_winkx(
            self.pointerToWinkx,
            self.pointerToListArrayWxscal,
            self.pointerToIwinkx,
            self.pointerToIdx,
            self.numberNodes,
            self.numberSlipperyWinklerForces,
            self.numberSlipperyWinklerEntries,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.slipperyWinklerInputFile,
            self.asciiOutputFile)
        # print "Just after read_winkx"
        # print "memorySize: %d" % self.memorySize
                             
        # except IOError, error:
            # print "Situation:", error
        # except ValueError, error:
            # print "Situation:", error
        # except ArithmeticError, error:
            # print "Situation:", error
        # except MemoryError, error:
            # print "Situation:", error
        # except Exception, error:
            # print "Exception in block read_winkx!", error
        
            
        # Localize global equation numbers in element index arrays.
        
        self.pointerToLm = lithomop3d.allocateInt(
            self.numberDegreesFreedom*self.connectivitySize)
	self.memorySize += self.numberDegreesFreedom* \
	    self.connectivitySize*self.intSize
        self.pointerToLmx = lithomop3d.allocateInt(
            self.numberDegreesFreedom*self.connectivitySize)
	self.memorySize += self.numberDegreesFreedom* \
	    self.connectivitySize*self.intSize
        # For now, try allocating everything whether it is used or not.
        #        if self.numberSplitNodeEntries != 0:
        #            self.pointerToLmf = lithomop3d.allocateInt(
        #                self.numberElementNodes*self.numberElements)
        self.pointerToLmf = lithomop3d.allocateInt(
            self.connectivitySize)
	self.memorySize += self.connectivitySize*self.intSize

        lithomop3d.local(
            self.pointerToId,
            self.numberNodes,
            self.pointerToIen,
            self.pointerToLm,
            self.pointerToInfiel,
            self.connectivitySize,
            self.numberElements,
            self.pointerToElementTypeInfo)
        # print "Just after local"
	# print "memorySize: %d" % self.memorySize

        lithomop3d.localf(
            self.pointerToIen,
            self.pointerToLmf,
            self.pointerToInfiel,
            self.connectivitySize,
            self.numberElements,
            self.pointerToElementTypeInfo,
            self.pointerToNfault,
            self.numberSplitNodeEntries)
        # print "Just after localf"

        lithomop3d.localx(
            self.pointerToIdx,
            self.numberNodes,
            self.pointerToIen,
            self.pointerToLmx,
            self.pointerToInfiel,
            self.connectivitySize,
            self.numberElements,
            self.pointerToElementTypeInfo,
            self.pointerToNslip,
            self.numberSlipperyNodeEntries)
        # print "Just after localx"
        self.pointerToNslip = None
	self.memorySize -= self.numberSlipDimensions* \
	    self.numberSlipperyNodeEntries*self.intSize

        # Allocate and populate sparse matrix arrays.  Some of these are
        # temporary and are then deleted after use.
        self.workingArraySize = lithomop3d.cmp_stiffsz(
            self.numberGlobalEquations,
            self.pointerToLm,
            self.pointerToLmx,
            self.pointerToInfiel,
            self.connectivitySize,
            self.numberElements,
            self.pointerToElementTypeInfo,
            self.totalNumberSlipperyNodes)

        # print "workingArraySize: %i" % self.workingArraySize

        # self.workingArraySize = 5000*self.numberGlobalEquations       # This may need to be adjusted.
        # self.workingArraySize = 100*self.numberGlobalEquations       # This may need to be adjusted.
        self.pointerToIndx = lithomop3d.allocateInt(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.intSize
        self.pointerToLink = lithomop3d.allocateInt(
            self.workingArraySize)
	self.memorySize += self.workingArraySize*self.intSize
        self.pointerToNbrs = lithomop3d.allocateInt(
            self.workingArraySize)
	self.memorySize += self.workingArraySize*self.intSize

        # print "Just before lnklst"
	# print "memorySize: %d" % self.memorySize
	# try:
        self.stiffnessMatrixInfo = lithomop3d.lnklst(
            self.numberGlobalEquations,
            self.pointerToLm,
            self.pointerToLmx,
            self.pointerToInfiel,
            self.connectivitySize,
            self.numberElements,
            self.pointerToElementTypeInfo,
            self.pointerToIndx,
            self.pointerToLink,
            self.pointerToNbrs,
            self.workingArraySize,
            self.totalNumberSlipperyNodes)
        # print "Just after lnklst"
        self.stiffnessMatrixSize = self.stiffnessMatrixInfo[0]
        self.stiffnessOffDiagonalSize = self.stiffnessMatrixInfo[1]
	self.stiffnessTrueSize = self.stiffnessMatrixSize-1
        # print "workingArraySize: %i" % self.workingArraySize
        # print "stiffnessMatrixSize: %i" % self.stiffnessMatrixSize
        # print "stiffnessOffDiagonalSize: %i" % self.stiffnessOffDiagonalSize
        # except IOError, error:
            # print "Situation:", error
        # except ValueError, error:
            # print "Situation:", error
        # except ArithmeticError, error:
            # print "Situation:", error
        # except MemoryError, error:
            # print "Situation:", error
        # except Exception, error:
            # print "Exception in block lnklst!", error

        self.pointerToJa = lithomop3d.allocateInt(
            self.stiffnessMatrixSize)
	self.memorySize += self.stiffnessMatrixSize*self.intSize

        # print "Just before makemsr"
	# print "stiffnessMatrixSize: %i" % self.stiffnessMatrixSize
	# print "stiffnessOffDiagonalSize: %i" % self.stiffnessOffDiagonalSize
	# print "memorySize: %d" % self.memorySize
        # lithomop3d.makemsr(
            # self.pointerToJa,
            # self.pointerToIndx,
            # self.pointerToLink,
            # self.pointerToNbrs,
            # self.numberGlobalEquations,
            # self.stiffnessMatrixSize,
            # self.workingArraySize)

        self.stiffnessMatrixStats = lithomop3d.makemsr(
            self.pointerToJa,
            self.pointerToIndx,
            self.pointerToLink,
            self.pointerToNbrs,
            self.numberGlobalEquations,
            self.stiffnessMatrixSize,
            self.workingArraySize)
        # print "Just after makemsr"
        self.minimumNonzeroTermsPerRow = self.stiffnessMatrixStats[0]
        self.maximumNonzeroTermsPerRow = self.stiffnessMatrixStats[1]
        self.averageNonzeroTermsPerRow = float(self.stiffnessMatrixStats[2])

        self.pointerToIndx = None
        self.pointerToLink = None
        self.pointerToNbrs = None
	self.memorySize -= self.numberGlobalEquations*self.intSize
	self.memorySize -= self.workingArraySize*self.intSize
	self.memorySize -= self.workingArraySize*self.intSize

        # Output sparse matrix info

        lithomop3d.write_sparse_info(
            self.numberGlobalEquations,
            self.stiffnessMatrixSize,
            self.minimumNonzeroTermsPerRow,
            self.maximumNonzeroTermsPerRow,
            self.averageNonzeroTermsPerRow,
            self.asciiOutputInt,
            self.f77AsciiOutput,
            self.asciiOutputFile)
        # print "Just after write_sparse_info"
        
        # Allocate memory for all additional arrays

        # Sparse matrix arrays
        self.pointerToAlnz = lithomop3d.allocateDouble(
            self.stiffnessMatrixSize)
	self.memorySize += self.stiffnessMatrixSize*self.doubleSize
        self.pointerToPcg = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToZcg = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize

        # Force vectors
        self.pointerToB = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToBtot = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToBres = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToGvec1 = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToGvec2 = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToPvec = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize

        # Displacement arrays
        self.pointerToD = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDeld = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDprev = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToDcur = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize

        # Slippery node arrays
        self.pointerToDx = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDeldx = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDxcur = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize

        # Split node arrays
        self.pointerToDfault = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSplitNodeEntries)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberSplitNodeEntries*self.doubleSize
        self.pointerToTfault = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSplitNodeEntries)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberSplitNodeEntries*self.doubleSize

        # Local stiffness matrix arrays
        self.pointerToS = lithomop3d.allocateDouble(
            self.maxElementEquations*self.maxElementEquations)
	self.memorySize += self.maxElementEquations* \
	    self.maxElementEquations*self.doubleSize
        self.pointerToStemp = lithomop3d.allocateDouble(
            self.maxElementEquations*self.maxElementEquations)
	self.memorySize += self.maxElementEquations* \
	    self.maxElementEquations*self.doubleSize

        # Element arrays
        self.pointerToState = lithomop3d.allocateDouble(
            self.stateVariableDimension*self.stateSize)
	self.memorySize += self.stateVariableDimension* \
	    self.stateSize*self.doubleSize
        self.pointerToDstate = lithomop3d.allocateDouble(
            self.stateVariableDimension*self.stateSize)
	self.memorySize += self.stateVariableDimension* \
	    self.stateSize*self.doubleSize
        self.pointerToDmat = lithomop3d.allocateDouble(
            self.materialMatrixDimension*self.dmatSize)
	self.memorySize += self.materialMatrixDimension* \
	    self.dmatSize*self.doubleSize
        self.pointerToListArrayIddmat = lithomop3d.intListToArray( 
            self.listIddmat)
	self.memorySize += 36*self.intSize
	# print "memorySize: %d" % self.memorySize


        # Create arrays from lists that will be needed for the solution

        # ncodat array
        self.listNcodat = [
            self.analysisTypeInt,
            self.debuggingOutputInt]
        self.pointerToListArrayNcodat = lithomop3d.intListToArray(
            self.listNcodat)
	self.memorySize += 2*self.intSize
            
        # npar array
        self.listNpar = [
            self.numberElements,
            self.numberMaterials,
            self.numberTractionBc,
            self.numberSlipperyNodeEntries,
            self.numberSplitNodeEntries,
            self.prestressAutoComputeInt,
            self.stateSize,
            self.dmatSize,
            self.connectivitySize,
            self.numberDifferentialForceEntries,
            self.quadratureOrderInt]
        self.pointerToListArrayNpar = lithomop3d.intListToArray(
            self.listNpar)
	self.memorySize += 11*self.intSize

        # nprint array
        self.listNprint = [
            self.numberFullOutputs,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.ucdOutputInt]
        self.pointerToListArrayNprint = lithomop3d.intListToArray(
            self.listNprint)
	self.memorySize += 4*self.intSize

        # nsiter array
        self.listNsiter = [
            self.preconditionerTypeInt,
            self.maxPcgIterations]
        self.pointerToListArrayNsiter = lithomop3d.intListToArray(
            self.listNsiter)
	self.memorySize += 2*self.intSize

        # nsysdat array
        self.listNsysdat = [
            self.numberNodes,
            self.numberGlobalEquations,
            self.stiffnessMatrixSize,
            self.numberRotationEntries,
            self.numberPrestressEntries,
            self.totalNumberSlipperyNodes,
            self.totalNumberSplitNodes,
            self.propertyListSize,
            self.numberWinklerForces,
            self.numberSlipperyWinklerForces,
            self.autoRotateSlipperyNodesInt]
        self.pointerToListArrayNsysdat = lithomop3d.intListToArray(
            self.listNsysdat)
	self.memorySize += 11*self.intSize

        # nunits array
        self.listNunits = [
            self.f77StandardInput,
            self.f77StandardOutput,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.f77UcdOutput]
        self.pointerToListArrayNunits = lithomop3d.intListToArray(
            self.listNunits)
	self.memorySize += 6*self.intSize

        # nvisdat array
        self.listNvisdat = [
            self.numberCycles,
            self.numberTimeStepGroups,
            self.totalNumberTimeSteps,
            self.numberLoadHistories]
        self.pointerToListArrayNvisdat = lithomop3d.intListToArray(
            self.listNvisdat)
	self.memorySize += 4*self.intSize
            
        # rgiter array
        self.listRgiter = [
            self.stressTolerance,
            self.minimumStrainPerturbation,
            self.initialStrainPerturbation]
        self.pointerToListArrayRgiter = lithomop3d.doubleListToArray(
            self.listRgiter)
	self.memorySize += 3*self.doubleSize

        # rtimdat array
        self.currentTimeStepSize = 0.0
        self.currentAlfaParameter = 0.0
        self.listRtimdat = [
            self.currentTimeStepSize,
            self.currentAlfaParameter,
            self.prestressAutoComputePoisson]
        self.pointerToListArrayRtimdat = lithomop3d.doubleListToArray(
            self.listRtimdat)
	self.memorySize += 3*self.doubleSize

	print ""
	print ""
        print "Sparse matrix information:"
	print ""
        print "workingArraySize:          %i" % self.workingArraySize
        print "stiffnessMatrixSize:       %i" % self.stiffnessTrueSize
        print "stiffnessOffDiagonalSize:  %i" % self.stiffnessOffDiagonalSize
        print "minimumNonzeroTermsPerRow: %i" % self.minimumNonzeroTermsPerRow
        print "maximumNonzeroTermsPerRow: %i" % self.maximumNonzeroTermsPerRow
        print "averageNonzeroTermsPerRow: %g" % self.averageNonzeroTermsPerRow

	# print "memorySize: %d" % self.memorySize
	print ""
	print ""
        print "Hello from lm3dsetup.run (end)!"
        return


    def __init__(self):
        Component.__init__(self, "lm3dsetup", "setup")

	print ""
        print "Hello from lm3dsetup.__init__!"

        return



# version
# $Id: Lithomop3d_setup.py,v 1.10 2004/08/25 01:55:26 willic3 Exp $

# End of file 
