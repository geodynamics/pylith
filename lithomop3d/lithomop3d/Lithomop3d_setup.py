#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2005  All Rights Reserved
#
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
        lm3dscan.preinitialize()

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

        self.stressTolerance = lm3dscan.stressTolerance
        self.minimumStrainPerturbation = lm3dscan.minimumStrainPerturbation
        self.initialStrainPerturbation = lm3dscan.initialStrainPerturbation

        self.usePreviousDisplacementFlag = lm3dscan.usePreviousDisplacementFlag

        self.gravityX = lm3dscan.gravityX
        self.gravityY = lm3dscan.gravityY
        self.gravityZ = lm3dscan.gravityZ

        self.prestressAutoComputePoisson = lm3dscan.prestressAutoComputePoisson
        self.prestressAutoComputeYoungs = lm3dscan.prestressAutoComputeYoungs
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

                                                                   
        # Initialize and define some integer parameters based on string
        # or logical parameters in python

        self.quadratureOrderInt = 0
        if lm3dscan.quadratureOrder == "Full":
            self.quadratureOrderInt = 1
        elif quadratureOrder == "Reduced":
            self.quadratureOrderInt = 2
        elif quadratureOrder == "Selective":
            self.quadratureOrderInt = 3
        else:
            self.quadratureOrderInt = 1

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
        self.maxElementEquations = lm3dscan._maxElementEquations
        self.maxNumberVolumeElementFamilies = lm3dscan._maxNumberVolumeElementFamilies
        self.pointerToListArrayNumberElementNodesBase = lm3dscan._pointerToListArrayNumberElementNodesBase

        # Invariant parameters related to material model
        self.maxMaterialModels = lm3dscan._maxMaterialModels
        self.maxStateVariables = lm3dscan._maxStateVariables
        self.maxState0Variables = lm3dscan._maxState0Variables
        self.pointerToMaterialModelInfo = lm3dscan._pointerToMaterialModelInfo

        # Parameters derived from values in the inventory or the category 2 parameters.
        self.analysisTypeInt = lm3dscan._analysisTypeInt
        self.prestressAutoComputeInt = lm3dscan._prestressAutoComputeInt
        self.prestressAutoChangeElasticPropsInt = lm3dscan._prestressAutoChangeElasticPropsInt

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
        self.pointerToListArrayPropertyList = lm3dscan._pointerToListArrayPropertyList

        self.numberVolumeElements = lm3dscan._numberVolumeElements
        self.numberVolumeElementFamilies = lm3dscan._numberVolumeElementFamilies
        self.volumeElementType = lm3dscan._volumeElementType
        self.pointerToVolumeElementFamilyList = lm3dscan._pointerToVolumeElementFamilyList

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
        from ElementTypeDef import ElementTypeDef
        import lithomop3d

        # This function performs a lot of memory allocation, and also bundles several
        # scalar values into lists.  It is assumed that the scalar values contained in
        # the lists will no longer be accessible (or useful) from within python once
        # they are passed as arrays.  The contents of the lists should be maintained as
        # global variables, however.

	print ""
        print "Hello from lm3dsetup.run (begin)!"
        print "Reading problem definition and allocating storage:"

        lithomop3d.PetscInitialize()
        self.autoprestrStage, self.elasticStage, self.viscousStage, self.iterateEvent = lithomop3d.setupPETScLogging()

        eltype=ElementTypeDef()

        # Initialize parameters that are defined in this function

        self.numberGlobalEquations = 0

        self.totalNumberSplitNodes = 0
        self.totalNumberSlipperyNodes = 0

        self.externFlag = 0
        self.tractionFlag = 0
        self.gravityFlag = 0
        self.concForceFlag = 0
        self.numberConcForces = 0
        self.prestressFlag = 0
        
        self.pointerToBextern = None
        self.pointerToBtraction = None
        self.pointerToBgravity = None
        self.pointerToBconcForce = None
        self.pointerToBintern = None
        self.pointerToBresid = None
        self.pointerToDispVec = None
        self.pointerToDprev = None
        self.listNforce = [0, 0, 0, 0, 0, 0]
        self.pointerToListArrayNforce = None

        self.listGrav = [0.0, 0.0, 0.0]
        self.pointerToListArrayGrav = None

        self.pointerToX = None
        self.pointerToD = None
        self.pointerToDeld = None
        self.pointerToDcur = None
        self.pointerToId = None
        self.pointerToIwink = None
        self.pointerToWink = None
        self.listNsysdat = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.pointerToListArrayNsysdat = None
        self.pointerToListArrayIddmat = None

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
        self.pointerToState0 = None
        self.state0Size = 0
        self.pointerToDmat = None
        self.pointerToIen = None
        self.pointerToLm = None
        self.pointerToLmx = None
        self.pointerToLmf = None
        self.pointerToIvfamily = None
        self.listNpar = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.pointerToListArrayNpar = None
        self.elementSizeInfo = [0, 0, 0, 0]
        self.pointerToIvftmp = None

        self.pointerToIelno = None
        self.pointerToIside = None
        self.pointerToIhistry = None
        self.pointerToPres = None
        self.pointerToPdir = None

        self.pointerToMhist = None
        self.propertySize = 0

        self.elementTypeInfo = [0, 0, 0, 0]
        self.pointerToListArrayElementTypeInfo = None
        self.pointerToSh = None
        self.pointerToShj = None
        self.pointerToGauss = None
        self.numberVolumeElementNodes = 0
        self.numberVolumeElementGaussPoints = 0
        self.numberVolumeElementEquations = 0
	self.connectivitySize = 0

        self.pointerToHistry = None
        self.listRtimdat = [0.0, 0.0, 0.0]
        self.pointerToListArrayRtimdat = None
        self.listNtimdat = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.currentTimeStep = 0
        self.currentIterationsBetweenReform = 0
        self.currentStepsBetweenReform = 0
        self.currentLargeDeformationFlag = 0
        self.currentMaximumIterations = 0
        self.currentNumberTotalIterations = 0
        self.currentNumberReforms = 0
        self.currentNumberTotalPcgIterations = 0
	self.reformFlagInt = 0
        self.pointerToListArrayNtimdat = None
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
        
        self.pointerToSkew = None

        self.pointerToIprint = None
        self.listNcodat = [0, 0]
        self.pointerToListArrayNcodat = None
        self.listNunits = [0, 0, 0, 0, 0]
        self.pointerToListArrayNunits = None
        self.listNprint = [0, 0, 0]
        self.pointerToListArrayNprint = None
        self.pointerToIstatout = None
        self.pointerToNstatout = None


        # Arrays that can be deallocated after use.
        self.pointerToTimes = None
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

        self.listGrav = [
            self.gravityX.value,
            self.gravityY.value,
            self.gravityZ.value]
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

        # Set up global integration for volume elements.

        eltype.getdef(
            self.volumeElementType,
            self.quadratureOrderInt,
            self.numberSpaceDimensions,
            self.numberDegreesFreedom)

        self.elementTypeInfo = eltype.elementTypeInfo
        self.pointerToSh = eltype.pointerToSh
        self.pointerToShj = eltype.pointerToShj
        self.pointerToGauss = eltype.pointerToGauss
        self.numberVolumeElementNodes = eltype.numberVolumeElementNodes
        self.numberVolumeElementGaussPoints = eltype.numberVolumeElementGaussPoints
        self.numberVolumeElementEquations = eltype.numberVolumeElementEquations
	self.connectivitySize = self.numberVolumeElements*self.numberVolumeElementNodes
        self.pointerToListArrayElementTypeInfo = lithomop3d.intListToArray(
            self.elementTypeInfo)
	self.memorySize += 4*self.intSize

        # Node-based info (coordinates, displacement arrays, BC, and skew BC).

        self.pointerToX = lithomop3d.allocateDouble(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes* \
	    self.doubleSize
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

        self.numberConcForces = lithomop3d.read_bc(
            self.pointerToBond,
            self.displacementScaleFactor,
            self.velocityScaleFactor,
            self.forceScaleFactor,
            self.pointerToIbond,
            self.numberNodes,
            self.numberBcEntries,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.bcInputFile,
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


        # Output stress computation and subiteration parameters.
        lithomop3d.write_strscomp(
            self.stressTolerance.value,
            self.minimumStrainPerturbation,
            self.initialStrainPerturbation,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)
        
        lithomop3d.write_subiter(
            self.usePreviousDisplacementFlag,
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
            3*self.maxStateVariables)
	self.memorySize += 3*self.maxStateVariables*self.intSize
        self.pointerToNstatout = lithomop3d.allocateInt(3)
	self.memorySize += 3*self.intSize

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
            self.pointerToNstatout,
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
            

        # Allocate and read info on material properties, connectivities, and prestresses
        lithomop3d.write_element_info(
            self.numberVolumeElements,
            self.quadratureOrderInt,
            self.prestressAutoComputeInt,
            self.prestressAutoChangeElasticPropsInt,
            self.prestressAutoComputePoisson,
            self.prestressAutoComputeYoungs.value,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        self.pointerToIen = lithomop3d.allocateInt(
            self.numberVolumeElementNodes*self.numberVolumeElements)
	self.memorySize += self.numberVolumeElementNodes* \
                           self.numberVolumeElements*self.intSize
        self.pointerToIvfamily = lithomop3d.allocateInt(
            5*self.numberVolumeElementFamilies)
        self.memorySize += 5*self.numberVolumeElementFamilies*self.intSize

        self.pointerToIvftmp = lithomop3d.allocateInt(
            self.numberVolumeElementFamilies)
        self.memorySize += self.numberVolumeElementFamilies*self.intSize

        self.pointerToIndxiel = lithomop3d.allocateInt(
            self.numberVolumeElements)
	self.memorySize += self.numberVolumeElements*self.intSize

        if self.numberPrestressEntries != 0 or self.prestressAutoComputeInt != 0:
            self.prestressFlag = 1

        self.elementSizeInfo = lithomop3d.read_connect(
            self.pointerToMaterialModelInfo,
            self.pointerToVolumeElementFamilyList,
            self.numberVolumeElementNodes,
            self.numberVolumeElementGaussPoints,
            self.pointerToIen,
            self.pointerToIvfamily,
            self.pointerToIvftmp,
            self.pointerToIndxiel,
            self.maxNumberVolumeElementFamilies,
            self.numberVolumeElementFamilies,
            self.prestressFlag,
            self.numberVolumeElements,
            self.numberNodes,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.connectivityInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)

        self.stateSize = self.elementSizeInfo[0]
        self.state0Size = self.elementSizeInfo[1]
        self.propertySize = self.elementSizeInfo[2]

        self.pointerToMhist = lithomop3d.allocateInt(
            self.propertySize)
	self.memorySize += self.propertySize*self.intSize

        lithomop3d.write_props(
            self.pointerToListArrayPropertyList,
            self.pointerToListArrayGrav,
            self.pointerToIvfamily,
            self.pointerToMaterialModelInfo,
            self.numberVolumeElementFamilies,
            self.propertySize,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputFile,
            self.plotOutputFile)

        lithomop3d.read_mathist(
            self.pointerToMhist,
            self.pointerToIvfamily,
            self.numberVolumeElementFamilies,
            self.propertySize,
            self.numberLoadHistories,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.materialHistoryInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)

        # write mesh info to UCD file, if requested
        if self.ucdOutputInt == 1:
            lithomop3d.write_ucd_mesh(
                self.pointerToX,
                self.numberNodes,
                self.pointerToIen,
                self.pointerToIvfamily,
                self.numberVolumeElements,
                self.numberVolumeElementFamilies,
                self.pointerToSh,
                self.numberVolumeElementNodes,
                self.numberVolumeElementGaussPoints,
                self.volumeElementType,
                self.pointerToIstatout,
                self.pointerToNstatout,
                self.f77UcdOutput,
                self.ucdOutputRoot)
                
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

        self.totalNumberSplitNodes = lithomop3d.read_split(
            self.pointerToFault,
            self.pointerToNfault,
            self.pointerToIndxiel,
            self.numberSplitNodeEntries,
            self.numberNodes,
            self.numberVolumeElements,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.splitNodeInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)

        self.totalNumberSlipperyNodes = lithomop3d.read_slip(
            self.pointerToNslip,
            self.pointerToIndxiel,
            self.numberSlipperyNodeEntries,
            self.numberNodes,
            self.autoRotateSlipperyNodesInt,
            self.numberVolumeElements,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.slipperyNodeInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)

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

        self.pointerToIndxiel = None
	self.memorySize -= self.numberVolumeElements*self.intSize


        # Create Idftn array for split nodes.  At present, this array is only used within
        # the routine itself to output split node info to a plot file.  It can probably be
        # eliminated in the near future.
        
        self.pointerToIdftn = lithomop3d.allocateInt(
            self.totalNumberSplitNodes)
	self.memorySize += self.totalNumberSplitNodes*self.intSize

        lithomop3d.id_split(
            self.pointerToNfault,
            self.pointerToIdftn,
            self.numberNodes,
            self.numberSplitNodeEntries,
            self.totalNumberSplitNodes,
            self.f77PlotOutput,
            self.plotOutputInt,
            self.plotOutputFile)

        # Deallocating this array since it does not appear to be used elsewhere.
        # It may be useful in the future, though.
        self.pointerToIdftn = None
	self.memorySize -= self.totalNumberSplitNodes*self.intSize

        # Determine global equations and store equation numbers in Id and Idx.
        self.pointerToId = lithomop3d.allocateInt(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes* \
	    self.intSize
        self.pointerToIdx = lithomop3d.allocateInt(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes*self.intSize
        self.pointerToIdslp = lithomop3d.allocateInt(
            self.numberNodes)
	self.memorySize += self.numberNodes*self.intSize

        self.numberGlobalEquations = lithomop3d.create_id(
            self.pointerToId,
            self.pointerToIdx,
            self.pointerToIbond,
            self.pointerToNslip,
            self.pointerToIdslp,
            self.numberSlipperyNodeEntries,
            self.numberNodes,
            self.totalNumberSlipperyNodes)

        self.pointerToIpslp = lithomop3d.allocateInt(
            self.numberSlipNeighbors*self.totalNumberSlipperyNodes)
        self.memorySize += self.numberSlipNeighbors* \
                           self.totalNumberSlipperyNodes*self.intSize


        # If there are slippery nodes and the auto-rotation option is selected, find
        # neighboring nodes on the fault so that a best-fit plane can be determined at
        # each node.
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

            self.pointerToXtmp = None
            self.pointerToItmp = None
            self.pointerToItmp1 = None
            self.pointerToItmp2 = None
	    self.memorySize -= self.totalNumberSlipperyNodes*self.doubleSize
	    self.memorySize -= self.totalNumberSlipperyNodes*self.intSize
	    self.memorySize -= self.totalNumberSlipperyNodes*self.intSize
	    self.memorySize -= self.totalNumberSlipperyNodes*self.intSize
                         
        # Read Winkler forces and slippery Winkler forces.  Also, equation info
        # is stored in Iwink and Iwinkx.
        # All input is finished after this section.
        self.pointerToIwink = lithomop3d.allocateInt(
            2*self.numberWinklerForces)
	self.memorySize += 2*self.numberWinklerForces*self.intSize
        self.pointerToWink = lithomop3d.allocateDouble(
            self.numberWinklerForces)
	self.memorySize += self.numberWinklerForces*self.doubleSize
        self.pointerToWinkx = lithomop3d.allocateDouble(
            self.numberSlipperyWinklerForces)
	self.memorySize += self.numberSlipperyWinklerForces*self.doubleSize
        self.pointerToIwinkx = lithomop3d.allocateInt(
            2*self.numberSlipperyWinklerForces)
	self.memorySize += 2*self.numberSlipperyWinklerForces*self.intSize

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
            self.numberVolumeElements,
            self.numberVolumeElementNodes)

        lithomop3d.localf(
            self.pointerToIen,
            self.pointerToLmf,
            self.numberVolumeElements,
            self.pointerToNfault,
            self.numberSplitNodeEntries,
            self.numberVolumeElementNodes)

        lithomop3d.localx(
            self.pointerToIdx,
            self.numberNodes,
            self.pointerToIen,
            self.pointerToLmx,
            self.numberVolumeElements,
            self.pointerToNslip,
            self.numberSlipperyNodeEntries,
            self.numberVolumeElementNodes)
        self.pointerToNslip = None
	self.memorySize -= self.numberSlipDimensions* \
	    self.numberSlipperyNodeEntries*self.intSize

        # Allocate and populate sparse matrix arrays.  Some of these are
        # temporary and are then deleted after use.
        self.workingArraySize = lithomop3d.cmp_stiffsz(
            self.numberGlobalEquations,
            self.pointerToLm,
            self.pointerToLmx,
            self.numberVolumeElements,
            self.totalNumberSlipperyNodes,
            self.numberVolumeElementNodes)


        self.pointerToIndx = lithomop3d.allocateInt(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.intSize
        self.pointerToLink = lithomop3d.allocateInt(
            self.workingArraySize)
	self.memorySize += self.workingArraySize*self.intSize
        self.pointerToNbrs = lithomop3d.allocateInt(
            self.workingArraySize)
	self.memorySize += self.workingArraySize*self.intSize

        self.stiffnessMatrixInfo = lithomop3d.lnklst(
            self.numberGlobalEquations,
            self.pointerToLm,
            self.pointerToLmx,
            self.numberVolumeElements,
            self.numberVolumeElementNodes,
            self.numberVolumeElementEquations,
            self.pointerToIndx,
            self.pointerToLink,
            self.pointerToNbrs,
            self.workingArraySize,
            self.totalNumberSlipperyNodes)

        self.stiffnessMatrixSize = self.stiffnessMatrixInfo[0]
        self.stiffnessOffDiagonalSize = self.stiffnessMatrixInfo[1]
	self.stiffnessTrueSize = self.stiffnessMatrixSize-1

        self.A = lithomop3d.createPETScMat(self.numberGlobalEquations)
	self.memorySize += self.stiffnessMatrixSize*self.intSize

        self.stiffnessMatrixStats = lithomop3d.makemsr(
            self.A,
            self.pointerToIndx,
            self.pointerToLink,
            self.pointerToNbrs,
            self.numberGlobalEquations,
            self.stiffnessMatrixSize,
            self.workingArraySize)

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
        
        # Allocate memory for all additional arrays

        # Sparse matrix arrays
	self.memorySize += self.stiffnessMatrixSize*self.doubleSize

        # Force vectors
        if self.numberTractionBc != 0:
            self.tractionFlag = 1
        if self.gravityX.value != 0.0 or self.gravityY.value != 0.0 or self.gravityZ.value != 0.0:
            self.gravityFlag = 1
        if self.numberConcForces != 0 or self.numberDifferentialForceEntries != 0:
            self.concForceFlag = 1
        if self.tractionFlag != 0 or self.gravityFlag != 0 or self.concForceFlag != 0:
            self.externFlag = 1

        self.pointerToBextern = lithomop3d.allocateDouble(
            self.externFlag*self.numberGlobalEquations)
	self.memorySize += self.externFlag*self.numberGlobalEquations*self.doubleSize
        self.pointerToBtraction = lithomop3d.allocateDouble(
            self.tractionFlag*self.numberGlobalEquations)
	self.memorySize += self.tractionFlag*self.numberGlobalEquations*self.doubleSize
        self.pointerToBgravity = lithomop3d.allocateDouble(
            self.gravityFlag*self.numberGlobalEquations)
	self.memorySize += self.gravityFlag*self.numberGlobalEquations*self.doubleSize
        self.pointerToBconcForce = lithomop3d.allocateDouble(
            self.concForceFlag*self.numberGlobalEquations)
	self.memorySize += self.concForceFlag*self.numberGlobalEquations*self.doubleSize
        self.pointerToBintern = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToBresid = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToDispVec = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToDprev = lithomop3d.allocateDouble(
            self.usePreviousDisplacementFlag*self.numberGlobalEquations)
	self.memorySize += self.usePreviousDisplacementFlag*self.numberGlobalEquations*self.doubleSize
            
        # Displacement arrays
        self.pointerToD = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDeld = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
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
            self.stateSize)
	self.memorySize += self.stateSize*self.doubleSize
        self.pointerToDstate = lithomop3d.allocateDouble(
            self.stateSize)
	self.memorySize += self.stateSize*self.doubleSize
        self.pointerToDmat = lithomop3d.allocateDouble(
            self.materialMatrixDimension*
            self.numberVolumeElementGaussPoints*
            self.numberVolumeElements)
	self.memorySize += self.materialMatrixDimension* \
	    self.numberVolumeElementGaussPoints* \
            self.numberVolumeElements*self.doubleSize
        self.pointerToListArrayIddmat = lithomop3d.intListToArray( 
            self.listIddmat)
	self.memorySize += 36*self.intSize
        self.pointerToState0 = lithomop3d.allocateDouble(
            self.state0Size)
        self.memorySize += self.state0Size*self.doubleSize

        # Create arrays from lists that will be needed for the solution

        # nforce array
        self.listNforce = [
            self.externFlag,
            self.tractionFlag,
            self.gravityFlag,
            self.concForceFlag,
            self.prestressFlag,
            self.usePreviousDisplacementFlag]
        self.pointerToListArrayNforce = lithomop3d.intListToArray(
            self.listNforce)
	self.memorySize += 6*self.intSize
           
        # ncodat array
        self.listNcodat = [
            self.analysisTypeInt,
            self.debuggingOutputInt]
        self.pointerToListArrayNcodat = lithomop3d.intListToArray(
            self.listNcodat)
	self.memorySize += 2*self.intSize
            
        # npar array
        self.listNpar = [
            self.numberVolumeElements,
            self.numberMaterials,
            self.numberTractionBc,
            self.numberSlipperyNodeEntries,
            self.numberSplitNodeEntries,
            self.prestressAutoComputeInt,
            self.prestressAutoChangeElasticPropsInt,
            self.stateSize,
            self.state0Size,
            self.numberVolumeElementFamilies,
            self.numberDifferentialForceEntries,
            self.quadratureOrderInt]
        self.pointerToListArrayNpar = lithomop3d.intListToArray(
            self.listNpar)
	self.memorySize += 12*self.intSize

        # nprint array
        self.listNprint = [
            self.numberFullOutputs,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.ucdOutputInt]
        self.pointerToListArrayNprint = lithomop3d.intListToArray(
            self.listNprint)
	self.memorySize += 4*self.intSize

        # nsysdat array
        self.listNsysdat = [
            self.numberNodes,
            self.numberGlobalEquations,
            self.stiffnessMatrixSize,
            self.numberRotationEntries,
            self.numberPrestressEntries,
            self.totalNumberSlipperyNodes,
            self.totalNumberSplitNodes,
            self.propertySize,
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
            self.stressTolerance.value,
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
            self.prestressAutoComputePoisson,
            self.prestressAutoComputeYoungs.value]
        self.pointerToListArrayRtimdat = lithomop3d.doubleListToArray(
            self.listRtimdat)
	self.memorySize += 4*self.doubleSize

        # ntimdat array
        self.listNtimdat = [
            self.currentTimeStep,
            self.currentIterationsBetweenReform,
            self.currentStepsBetweenReform,
            self.currentLargeDeformationFlag,
            self.currentMaximumIterations,
            self.currentNumberTotalIterations,
            self.currentNumberReforms,
            self.currentNumberTotalPcgIterations,
            self.reformFlagInt]
        self.pointerToListArrayNtimdat = lithomop3d.intListToArray(
            self.listNtimdat)
        self.memorySize += 9*self.intSize

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
# $Id: Lithomop3d_setup.py,v 1.24 2005/04/08 17:38:16 willic3 Exp $

# End of file 
