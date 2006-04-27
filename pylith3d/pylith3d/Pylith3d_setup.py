#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams
#  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
#  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
#  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
# The code here should be executed after all initializations in Pylith3d_scan.py
# have been performed, including reading of the keyword=value file.
#


from pyre.components.Component import Component


class Pylith3d_setup(Component):

    def initialize(self, scanner):

        pl3dscan = scanner
        pl3dscan.preinitialize()

	print ""
        print "Hello from pl3dsetup.initialize (begin)!"
        print "Importing values from scanning phase:"

        # Parameters needed from Pylith3d_scan.py.
        # These parameters, which may be modified in the run function, should all be
        # available to Pylith3d_run.py.  Inventory items are not altered here and
        # may be accessed from Pylith3d_scan.py.

        self.winklerScaleX = pl3dscan.winklerScaleX
        self.winklerScaleY = pl3dscan.winklerScaleY
        self.winklerScaleZ = pl3dscan.winklerScaleZ

        self.stressTolerance = pl3dscan.stressTolerance
        self.minimumStrainPerturbation = pl3dscan.minimumStrainPerturbation
        self.initialStrainPerturbation = pl3dscan.initialStrainPerturbation

        self.usePreviousDisplacementFlag = pl3dscan.usePreviousDisplacementFlag

        self.gravityX = pl3dscan.gravityX
        self.gravityY = pl3dscan.gravityY
        self.gravityZ = pl3dscan.gravityZ

        self.prestressAutoComputePoisson = pl3dscan.prestressAutoComputePoisson
        self.prestressAutoComputeYoungs = pl3dscan.prestressAutoComputeYoungs
        self.prestressScaleXx = pl3dscan.prestressScaleXx
        self.prestressScaleYy = pl3dscan.prestressScaleYy
        self.prestressScaleZz = pl3dscan.prestressScaleZz
        self.prestressScaleXy = pl3dscan.prestressScaleXy
        self.prestressScaleXz = pl3dscan.prestressScaleXz
        self.prestressScaleYz = pl3dscan.prestressScaleYz

        self.winklerSlipScaleX = pl3dscan.winklerSlipScaleX
        self.winklerSlipScaleY = pl3dscan.winklerSlipScaleY
        self.winklerSlipScaleZ = pl3dscan.winklerSlipScaleZ

        self.f77StandardInput = pl3dscan.f77StandardInput
        self.f77StandardOutput = pl3dscan.f77StandardOutput
        self.f77FileInput = pl3dscan.f77FileInput
        self.f77AsciiOutput = pl3dscan.f77AsciiOutput
        self.f77PlotOutput = pl3dscan.f77PlotOutput
        self.f77UcdOutput = pl3dscan.f77UcdOutput

                                                                   
        # Initialize and define some integer parameters based on string
        # or logical parameters in python

        self.quadratureOrderInt = 0
        if pl3dscan.quadratureOrder == "Full":
            self.quadratureOrderInt = 1
        elif quadratureOrder == "Reduced":
            self.quadratureOrderInt = 2
        elif quadratureOrder == "Selective":
            self.quadratureOrderInt = 3
        else:
            self.quadratureOrderInt = 1

        self.asciiOutputInt = 0
        if pl3dscan.inventory.asciiOutput == "none":
            self.asciiOutputInt = 0
        elif pl3dscan.inventory.asciiOutput == "echo":
            self.asciiOutputInt = 1
        else:
            self.asciiOutputInt = 2
            
        self.plotOutputInt = 0
        if pl3dscan.inventory.plotOutput == "none":
            self.plotOutputInt = 0
        elif pl3dscan.inventory.plotOutput == "ascii":
            self.plotOutputInt = 1
        else:
            self.plotOutputInt = 2

        binIOError = None
        import pylith3d
        try:
            pylith3d.try_binio(self.f77UcdOutput)
        except RuntimeError, binIOError:
            self.ucdOutputInt = 1
        else:
            self.ucdOutputInt = 2
        if pl3dscan.inventory.ucdOutput == "none":
            self.ucdOutputInt = 0
        elif pl3dscan.inventory.ucdOutput == "ascii":
            self.ucdOutputInt = 1
        elif pl3dscan.inventory.ucdOutput == "binary":
            if binIOError is None:
                self.ucdOutputInt = 2
            else:
                import journal
                warning = journal.warning("pylith3d")
                warning.line("Forcing 'pl3dscan.ucdOutput' to 'ascii'.")
                warning.line("Binary UCD output not supported for this Fortran compiler.")
                warning.log(binIOError)
            
        self.debuggingOutputInt = 0
        if pl3dscan.inventory.debuggingOutput:
            self.debuggingOutputInt = 1
        else:
            self.debuggingOutputInt = 0

        self.autoRotateSlipperyNodesInt = 0
        if pl3dscan.inventory.autoRotateSlipperyNodes:
            self.autoRotateSlipperyNodesInt = 2
        else:
            self.autoRotateSlipperyNodesInt = 1

        # Get some parameters from the inventory list.
        self.title = pl3dscan.inventory.title
        self.numberCycles = pl3dscan.inventory.numberCycles

        # Category 2 parameters needed from pl3dscan._init.  These parameters are
        # not meant to be altered by the user.
        
	# Approximate memory allocation info:
	self.memorySize = pl3dscan._memorySize
	self.intSize = pl3dscan._intSize
	self.doubleSize = pl3dscan._doubleSize

        # First get all filenames
        self.asciiOutputFile = pl3dscan._asciiOutputFile
        self.plotOutputFile = pl3dscan._plotOutputFile
        self.ucdOutputRoot = pl3dscan._ucdOutputRoot
        self.coordinateInputFile = pl3dscan._coordinateInputFile
        self.bcInputFile = pl3dscan._bcInputFile
        self.winklerInputFile = pl3dscan._winklerInputFile
        self.rotationInputFile = pl3dscan._rotationInputFile
        self.timeStepInputFile = pl3dscan._timeStepInputFile
        self.fullOutputInputFile = pl3dscan._fullOutputInputFile
        self.stateVariableInputFile = pl3dscan._stateVariableInputFile
        self.loadHistoryInputFile = pl3dscan._loadHistoryInputFile
        self.materialPropertiesInputFile = pl3dscan._materialPropertiesInputFile
        self.materialHistoryInputFile = pl3dscan._materialHistoryInputFile
        self.connectivityInputFile = pl3dscan._connectivityInputFile
        self.prestressInputFile = pl3dscan._prestressInputFile
        self.tractionInputFile = pl3dscan._tractionInputFile
        self.splitNodeInputFile = pl3dscan._splitNodeInputFile
        self.slipperyNodeInputFile = pl3dscan._slipperyNodeInputFile
        self.differentialForceInputFile = pl3dscan._differentialForceInputFile
        self.slipperyWinklerInputFile = pl3dscan._slipperyWinklerInputFile

        # Parameters that are invariant for this geometry type
        self.geometryType = pl3dscan._geometryType
        self.geometryTypeInt = pl3dscan._geometryTypeInt
        self.numberSpaceDimensions = pl3dscan._numberSpaceDimensions
        self.numberDegreesFreedom = pl3dscan._numberDegreesFreedom
        self.stateVariableDimension = pl3dscan._stateVariableDimension
        self.materialMatrixDimension = pl3dscan._materialMatrixDimension
        self.numberSkewDimensions = pl3dscan._numberSkewDimensions
        self.numberSlipDimensions = pl3dscan._numberSlipDimensions
        self.numberSlipNeighbors = pl3dscan._numberSlipNeighbors
        self.numberTractionDirections = pl3dscan._numberTractionDirections
        self.listIddmat = pl3dscan._listIddmat

        # Invariant parameters related to element type
        self.maxElementEquations = pl3dscan._maxElementEquations
        self.maxNumberVolumeElementFamilies = pl3dscan._maxNumberVolumeElementFamilies
        self.pointerToListArrayNumberElementNodesBase = pl3dscan._pointerToListArrayNumberElementNodesBase

        # Invariant parameters related to material model
        self.maxMaterialModels = pl3dscan._maxMaterialModels
        self.maxStateVariables = pl3dscan._maxStateVariables
        self.maxState0Variables = pl3dscan._maxState0Variables
        self.pointerToMaterialModelInfo = pl3dscan._pointerToMaterialModelInfo

        # Parameters derived from values in the inventory or the category 2 parameters.
        self.analysisTypeInt = pl3dscan._analysisTypeInt
        self.prestressAutoComputeInt = pl3dscan._prestressAutoComputeInt
        self.prestressAutoChangeElasticPropsInt = pl3dscan._prestressAutoChangeElasticPropsInt

        # Parameters derived from the number of entries in a file
        self.numberNodes = pl3dscan._numberNodes
        self.coordinateScaleFactor = pl3dscan._coordinateScaleFactor

        self.numberBcEntries = pl3dscan._numberBcEntries
        self.displacementScaleFactor = pl3dscan._displacementScaleFactor
        self.velocityScaleFactor = pl3dscan._velocityScaleFactor
        self.forceScaleFactor = pl3dscan._forceScaleFactor

        self.numberWinklerEntries = pl3dscan._numberWinklerEntries
        self.numberWinklerForces = pl3dscan._numberWinklerForces

        self.numberRotationEntries = pl3dscan._numberRotationEntries
        self.rotationScaleFactor = pl3dscan._rotationScaleFactor

        self.numberTimeStepGroups = pl3dscan._numberTimeStepGroups
        self.totalNumberTimeSteps = pl3dscan._totalNumberTimeSteps
        self.timeScaleFactor = pl3dscan._timeScaleFactor

        self.numberFullOutputs = pl3dscan._numberFullOutputs

        self.numberLoadHistories = pl3dscan._numberLoadHistories

        self.numberMaterials = pl3dscan._numberMaterials
        self.pointerToListArrayPropertyList = pl3dscan._pointerToListArrayPropertyList

        self.numberVolumeElements = pl3dscan._numberVolumeElements
        self.numberVolumeElementFamilies = pl3dscan._numberVolumeElementFamilies
        self.volumeElementType = pl3dscan._volumeElementType
        self.pointerToVolumeElementFamilyList = pl3dscan._pointerToVolumeElementFamilyList

        self.numberPrestressEntries = pl3dscan._numberPrestressEntries

        self.numberTractionBc = pl3dscan._numberTractionBc
        self.tractionBcScaleFactor = pl3dscan._tractionBcScaleFactor

        self.numberSplitNodeEntries = pl3dscan._numberSplitNodeEntries

        self.numberSlipperyNodeEntries = pl3dscan._numberSlipperyNodeEntries
        self.numberDifferentialForceEntries = pl3dscan._numberDifferentialForceEntries
        self.numberSlipperyWinklerEntries = pl3dscan._numberSlipperyWinklerEntries
        self.numberSlipperyWinklerForces = pl3dscan._numberSlipperyWinklerForces

	print ""
        print "Hello from pl3dsetup.initialize (end)!"

        return

    def read(self):

        # This function reads all input and performs some memory allocation.

        from ElementTypeDef import ElementTypeDef
        import pylith3d

	print ""
        print "Hello from pl3dsetup.read (begin)!"
        print "Reading problem definition and allocating necessary storage:"


        eltype=ElementTypeDef()

        # Initialize variables that are defined in this function.

        # Number of split/slippery nodes
        self.totalNumberSplitNodes = 0
        self.totalNumberSlipperyNodes = 0

        # Force vector flags
        self.externFlag = 0
        self.tractionFlag = 0
        self.gravityFlag = 0
        self.concForceFlag = 0
        self.numberConcForces = 0
        self.prestressFlag = 0
	self.winklerFlag = 0
	self.slipperyWinklerFlag = 0

        # Nodal arrays and equation numbers
        self.pointerToX = None
        self.pointerToIwinkdef = None
        self.pointerToIwinkid = None
        self.pointerToWinkdef = None

        # Local coordinate rotations
        self.pointerToSkew = None

        # Nodal boundary conditions
        self.pointerToIbond = None
        self.pointerToBond = None

	# Time step information
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

        # Split node arrays
        self.pointerToFault = None
        self.pointerToNfault = None

        # Slippery node arrays
        self.pointerToDiforc = None
        self.pointerToIwinkxdef = None
        self.pointerToIwinkxid = None
        self.pointerToWinkxdef = None
        self.pointerToIdhist = None

        # Element information
        self.pointerToIen = None
        self.pointerToMat = None

        # Element type information
        self.elementTypeInfo = [0, 0, 0, 0]
        self.pointerToListArrayElementTypeInfo = None
        self.pointerToSh = None
        self.pointerToShj = None
        self.pointerToGauss = None
        self.numberVolumeElementNodes = 0
        self.numberVolumeElementGaussPoints = 0
        self.numberVolumeElementEquations = 0
	self.connectivitySize = 0

        # Traction BC
        self.pointerToIelno = None
        self.pointerToIside = None
        self.pointerToIhistry = None
        self.pointerToPres = None
        self.pointerToPdir = None

        # Time histories
        self.pointerToHistry = None

        # Output information
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
        # Note that array Nslip is also required in functions sortmesh and sparsesetup
        # before it can be deallocated.  Also, array Times is needed for output, if
        # requested.
        self.pointerToTimes = None
        self.pointerToNslip = None
        self.pointerToXtmp = None
        self.pointerToItmp = None
        self.pointerToItmp1 = None
        self.pointerToItmp2 = None

        # Lists that are used as arrays in the input routines below
        self.listWscal = [0.0, 0.0, 0.0]
        self.pointerToListArrayWscal = None
        self.listPrscal = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.pointerToListArrayPrscal = None
        self.listWxscal = [0.0, 0.0, 0.0]
        self.pointerToListArrayWxscal = None


        # Make lists that are used as arrays in the f77 function calls below.
        self.listWscal = [
            self.winklerScaleX,
            self.winklerScaleY,
            self.winklerScaleZ]
        self.pointerToListArrayWscal = pylith3d.doubleListToArray(
            self.listWscal)
	self.memorySize += 3*self.doubleSize

        self.listPrscal = [
            self.prestressScaleXx,
            self.prestressScaleYy,
            self.prestressScaleZz,
            self.prestressScaleXy,
            self.prestressScaleXz,
            self.prestressScaleYz]
        self.pointerToListArrayPrscal = pylith3d.doubleListToArray(
            self.listPrscal)
	self.memorySize += 6*self.doubleSize
                                  
        self.listWxscal = [
            self.winklerSlipScaleX,
            self.winklerSlipScaleY,
            self.winklerSlipScaleZ]
        self.pointerToListArrayWxscal = pylith3d.doubleListToArray(
            self.listWxscal)
	self.memorySize += 3*self.doubleSize

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
        self.pointerToListArrayElementTypeInfo = pylith3d.intListToArray(
            self.elementTypeInfo)
	self.memorySize += 4*self.intSize

        # Node-based info (coordinates, displacement arrays, BC, and skew BC).
        self.pointerToX = pylith3d.allocateDouble(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes* \
	    self.doubleSize
        self.pointerToIbond = pylith3d.allocateInt(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes* \
	    self.intSize
        self.pointerToBond = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes* \
	    self.doubleSize
        self.pointerToSkew = pylith3d.allocateDouble(
            self.numberSkewDimensions*self.numberNodes)
	self.memorySize += self.numberSkewDimensions* \
	    self.numberNodes* \
	    self.doubleSize

        pylith3d.read_coords(
            self.pointerToX,
            self.coordinateScaleFactor,
            self.numberNodes,
            self.f77FileInput,
            self.coordinateInputFile)

        self.numberConcForces = pylith3d.read_bc(
            self.pointerToBond,
            self.displacementScaleFactor,
            self.velocityScaleFactor,
            self.forceScaleFactor,
            self.pointerToIbond,
            self.numberNodes,
            self.numberBcEntries,
            self.f77FileInput,
            self.bcInputFile)

        pylith3d.read_skew(
            self.pointerToSkew,
            self.rotationScaleFactor,
            self.numberRotationEntries,
            self.numberNodes,
            self.autoRotateSlipperyNodesInt,
            self.f77FileInput,
            self.rotationInputFile)

        # Allocate and read time step, time output, and load history info.
        self.pointerToHistry = pylith3d.allocateDouble(
            (self.totalNumberTimeSteps+1)*self.numberLoadHistories)
	self.memorySize += (self.totalNumberTimeSteps+1)* \
	    self.numberLoadHistories* \
	    self.doubleSize
        self.pointerToMaxstp = pylith3d.allocateInt(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.intSize
        self.pointerToDelt = pylith3d.allocateDouble(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.doubleSize
        self.pointerToAlfa = pylith3d.allocateDouble(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.doubleSize
        self.pointerToMaxit = pylith3d.allocateInt(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.intSize
        self.pointerToNtdinit = pylith3d.allocateInt(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.intSize
        self.pointerToLgdef = pylith3d.allocateInt(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.intSize
        self.pointerToUtol = pylith3d.allocateDouble(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.doubleSize
        self.pointerToFtol = pylith3d.allocateDouble(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.doubleSize
        self.pointerToEtol = pylith3d.allocateDouble(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.doubleSize
        self.pointerToItmax = pylith3d.allocateInt(
            self.numberTimeStepGroups)
	self.memorySize += self.numberTimeStepGroups*self.intSize
        self.pointerToIprint = pylith3d.allocateInt(
            self.numberFullOutputs)
	self.memorySize += self.numberFullOutputs*self.intSize
        self.pointerToTimes = pylith3d.allocateDouble(
            self.totalNumberTimeSteps+1)
	self.memorySize += (self.totalNumberTimeSteps+1)*self.doubleSize
        self.pointerToIstatout = pylith3d.allocateInt(
            3*self.maxStateVariables)
	self.memorySize += 3*self.maxStateVariables*self.intSize
        self.pointerToNstatout = pylith3d.allocateInt(3)
	self.memorySize += 3*self.intSize

        pylith3d.read_timdat(
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
            self.timeStepInputFile)

        pylith3d.read_fuldat(
            self.pointerToIprint,
            self.numberFullOutputs,
            self.analysisTypeInt,
            self.numberCycles,
            self.totalNumberTimeSteps,
            self.f77FileInput,
            self.fullOutputInputFile)

        pylith3d.read_stateout(
            self.pointerToIstatout,
            self.pointerToNstatout,
            self.f77FileInput,
            self.stateVariableInputFile)

        pylith3d.read_hist(
            self.pointerToHistry,
            self.pointerToTimes,
            self.numberLoadHistories,
            self.totalNumberTimeSteps,
            self.f77FileInput,
            self.loadHistoryInputFile)

        # Allocate and read info on connectivities and prestresses
        self.pointerToIen = pylith3d.allocateInt(
            self.numberVolumeElementNodes*self.numberVolumeElements)
	self.memorySize += self.numberVolumeElementNodes* \
                           self.numberVolumeElements*self.intSize
	self.pointerToMat = pylith3d.allocateInt(
	    self.numberVolumeElements)
        self.memorySize += self.numberVolumeElements*self.intSize
        if self.numberPrestressEntries != 0 or self.prestressAutoComputeInt != 0:
            self.prestressFlag = 1

        pylith3d.read_connect(
            self.pointerToIen,
            self.pointerToMat,
            self.numberVolumeElementNodes,
            self.numberVolumeElements,
            self.numberNodes,
	    self.numberVolumeElementFamilies,
            self.f77FileInput,
            self.connectivityInputFile)

        # pylith3d.read_prestr(
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

        # Read traction BC (unimplemented for now).
        self.pointerToIelno = pylith3d.allocateInt(
            self.numberTractionBc)
        self.pointerToIside = pylith3d.allocateInt(
            self.numberTractionBc)
        self.pointerToIhistry = pylith3d.allocateInt(
            self.numberTractionBc)
	#  Note that the following dimension is definitely wrong.
        self.pointerToPres = pylith3d.allocateDouble(
            self.numberTractionBc)
        self.pointerToPdir = pylith3d.allocateDouble(
            self.numberTractionDirections*self.numberTractionBc)

        # pylith3d.read_traction(
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

        # Read split node info
        self.pointerToNfault = pylith3d.allocateInt(
            3*self.numberSplitNodeEntries)
	self.memorySize += 3*self.numberSplitNodeEntries*self.intSize
        self.pointerToFault = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSplitNodeEntries)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberSplitNodeEntries*self.doubleSize

        self.totalNumberSplitNodes = pylith3d.read_split(
            self.pointerToFault,
            self.pointerToNfault,
            self.numberSplitNodeEntries,
            self.numberNodes,
            self.numberVolumeElements,
            self.f77FileInput,
            self.splitNodeInputFile)

        # Read slippery node info
        self.pointerToNslip = pylith3d.allocateInt(
            self.numberSlipDimensions*self.numberSlipperyNodeEntries)
	self.memorySize += self.numberSlipDimensions* \
	    self.numberSlipperyNodeEntries*self.intSize
        self.pointerToIdhist = pylith3d.allocateInt(
            self.numberNodes)
	self.memorySize += self.numberNodes*self.intSize
        self.pointerToDiforc = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize

        self.totalNumberSlipperyNodes = pylith3d.read_slip(
            self.pointerToNslip,
            self.numberSlipperyNodeEntries,
            self.numberNodes,
            self.autoRotateSlipperyNodesInt,
            self.f77FileInput,
            self.slipperyNodeInputFile)

        pylith3d.read_diff(
            self.pointerToDiforc,
            self.pointerToNslip,
            self.pointerToIdhist,
            self.numberSlipperyNodeEntries,
            self.numberDifferentialForceEntries,
            self.numberNodes,
            self.f77FileInput,
            self.differentialForceInputFile)
                         
        # Read Winkler forces and slippery Winkler forces.
        # All input is finished after this section.
        self.pointerToIwinkdef = pylith3d.allocateInt(
            self.numberDegreesFreedom*self.numberWinklerEntries)
	self.memorySize += self.numberDegreesFreedom* \
                           self.numberWinklerEntries*self.intSize
        self.pointerToIwinkid = pylith3d.allocateInt(
            self.numberWinklerEntries)
	self.memorySize += self.numberWinklerEntries*self.intSize
        self.pointerToWinkdef = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberWinklerEntries)
	self.memorySize += self.numberDegreesFreedom* \
                           self.numberWinklerEntries*self.doubleSize

        self.pointerToIwinkxdef = pylith3d.allocateInt(
            self.numberDegreesFreedom*self.numberSlipperyWinklerEntries)
	self.memorySize += self.numberDegreesFreedom* \
                           self.numberSlipperyWinklerEntries*self.intSize
        self.pointerToIwinkxid = pylith3d.allocateInt(
            self.numberSlipperyWinklerEntries)
	self.memorySize += self.numberSlipperyWinklerEntries*self.intSize
        self.pointerToWinkxdef = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSlipperyWinklerEntries)
	self.memorySize += self.numberDegreesFreedom* \
                           self.numberSlipperyWinklerEntries*self.doubleSize

        pylith3d.read_wink(
            self.pointerToWinkdef,
            self.pointerToListArrayWscal,
            self.pointerToIwinkdef,
            self.pointerToIwinkid,
            self.numberWinklerForces,
            self.numberWinklerEntries,
            self.f77FileInput,
            self.winklerInputFile)

        pylith3d.read_wink(
            self.pointerToWinkxdef,
            self.pointerToListArrayWxscal,
            self.pointerToIwinkxdef,
            self.pointerToIwinkxid,
            self.numberSlipperyWinklerForces,
            self.numberSlipperyWinklerEntries,
            self.f77FileInput,
            self.slipperyWinklerInputFile)

	print ""
        print "Hello from pl3dsetup.read (end)!"

        return

    def numberequations(self):

        # This functions numbers equations based on BC and slippery node info.

        import pylith3d

	print ""
        print "Hello from pl3dsetup.numberequations (begin)!"
        print "Numbering global equations:"

        # Initialize variables that are defined in this function.
        
        # Number of equations
        self.numberGlobalEquations = 0

        # Nodal equation numbers and Winkler restoring force info
        self.pointerToId = None
        self.pointerToIwink = None
        self.pointerToWink = None

        # Split node ID array.  This can be deallocated after meshwrite function has been called.
        self.pointerToIdftn = None

        # Slippery node equation numbers and Winkler restoring force info
        self.pointerToIdx = None
        self.pointerToIwinkx = None
        self.pointerToWinkx = None
        self.pointerToIdslp = None
        self.pointerToIpslp = None

        # Create Idftn array for split nodes.
        self.pointerToIdftn = pylith3d.allocateInt(
            self.totalNumberSplitNodes)
	self.memorySize += self.totalNumberSplitNodes*self.intSize

        pylith3d.id_split(
            self.pointerToNfault,
            self.pointerToIdftn,
            self.numberNodes,
            self.numberSplitNodeEntries,
            self.totalNumberSplitNodes)

        # Determine global equations and store equation numbers in Id and Idx.
        self.pointerToId = pylith3d.allocateInt(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes* \
	    self.intSize
        self.pointerToIdx = pylith3d.allocateInt(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes*self.intSize
        self.pointerToIdslp = pylith3d.allocateInt(
            self.numberNodes)
	self.memorySize += self.numberNodes*self.intSize

        self.numberGlobalEquations = pylith3d.create_id(
            self.pointerToId,
            self.pointerToIdx,
            self.pointerToIbond,
            self.pointerToNslip,
            self.pointerToIdslp,
            self.numberSlipperyNodeEntries,
            self.numberNodes,
            self.totalNumberSlipperyNodes)

        self.pointerToIpslp = pylith3d.allocateInt(
            self.numberSlipNeighbors*self.totalNumberSlipperyNodes)
        self.memorySize += self.numberSlipNeighbors* \
                           self.totalNumberSlipperyNodes*self.intSize

        # If there are slippery nodes and the auto-rotation option is selected, find
        # neighboring nodes on the fault so that a best-fit plane can be determined at
        # each node.  Deallocate temporary arrays after use.
        if self.totalNumberSlipperyNodes != 0 and self.autoRotateSlipperyNodesInt ==  2:

            self.pointerToXtmp = pylith3d.allocateDouble(
                self.totalNumberSlipperyNodes)
	    self.memorySize += self.totalNumberSlipperyNodes*self.doubleSize
            self.pointerToItmp = pylith3d.allocateInt(
                self.totalNumberSlipperyNodes)
	    self.memorySize += self.totalNumberSlipperyNodes*self.intSize
            self.pointerToItmp1 = pylith3d.allocateInt(
                self.totalNumberSlipperyNodes)
	    self.memorySize += self.totalNumberSlipperyNodes*self.intSize
            self.pointerToItmp2 = pylith3d.allocateInt(
                self.totalNumberSlipperyNodes)
	    self.memorySize += self.totalNumberSlipperyNodes*self.intSize

            pylith3d.nfind(
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

        # Assign appropriate equation numbers to Iwink array, and compact Wink
        # array to correspond to assigned BC.
        self.pointerToWink = pylith3d.allocateDouble(
            self.numberWinklerForces)
        self.memorySize += self.numberWinklerForces*self.doubleSize
        self.pointerToIwink = pylith3d.allocateInt(
            2*self.numberWinklerForces)
        self.memorySize += 2*self.numberWinklerForces*self.intSize

        pylith3d.assign_wink(
            self.pointerToWinkdef,
            self.pointerToWink,
            self.pointerToIwinkdef,
            self.pointerToIwinkid,
            self.pointerToIwink,
            self.pointerToId,
            self.numberNodes,
            self.numberWinklerForces,
            self.numberWinklerEntries)

        # Assign appropriate equation numbers to Iwinkx array, and compact Winkx
        # array to correspond to assigned BC.
        self.pointerToWinkx = pylith3d.allocateDouble(
            self.numberSlipperyWinklerForces)
        self.memorySize += self.numberSlipperyWinklerForces*self.doubleSize
        self.pointerToIwinkx = pylith3d.allocateInt(
            2*self.numberSlipperyWinklerForces)
        self.memorySize += 2*self.numberSlipperyWinklerForces*self.intSize

        pylith3d.assign_wink(
            self.pointerToWinkxdef,
            self.pointerToWinkx,
            self.pointerToIwinkxdef,
            self.pointerToIwinkxid,
            self.pointerToIwinkx,
            self.pointerToIdx,
            self.numberNodes,
            self.numberSlipperyWinklerForces,
            self.numberSlipperyWinklerEntries)

	print ""
        print "Hello from pl3dsetup.numberequations (end)!"
            
        return
            
        
    def sortmesh(self):

        # This function sorts elements into families and sorts all other items that are
        # affected by this.

        import pylith3d

	print ""
        print "Hello from pl3dsetup.sortmesh (begin)!"
        print "Renumbering elements, traction BC, split nodes, and slippery nodes:"

        # Initialize variables that are defined in this function.

        # Element arrays
        self.pointerToIens = None
        self.pointerToIvfamily = None
        self.elementSizeInfo = [ 0, 0, 0]
        self.pointerToIvftmp = None
        self.stateSize = 0
        self.state0Size = 0
        self.propertySize = 0

        # Sort elements into families.  The sorted elements are contained
        # in array Iens, and the index array for the original ordering is
        # Indxiel.  The original element node array (Ien) and the associated
        # material type array (Mat) may be deallocated after sorting.
        self.pointerToIens = pylith3d.allocateInt(
            self.numberVolumeElementNodes*self.numberVolumeElements)
	self.memorySize += self.numberVolumeElementNodes* \
                           self.numberVolumeElements*self.intSize
        self.pointerToIvfamily = pylith3d.allocateInt(
            5*self.numberVolumeElementFamilies)
        self.memorySize += 5*self.numberVolumeElementFamilies*self.intSize

        self.pointerToIvftmp = pylith3d.allocateInt(
            self.numberVolumeElementFamilies)
        self.memorySize += self.numberVolumeElementFamilies*self.intSize

        self.pointerToIndxiel = pylith3d.allocateInt(
            self.numberVolumeElements)
	self.memorySize += self.numberVolumeElements*self.intSize

	self.elementSizeInfo = pylith3d.sort_elements(
            self.pointerToIen,
            self.pointerToMat,
            self.pointerToMaterialModelInfo,
            self.pointerToVolumeElementFamilyList,
            self.pointerToIvfamily,
            self.pointerToIens,
            self.pointerToIvftmp,
            self.pointerToIndxiel,
            self.numberVolumeElementNodes,
            self.numberVolumeElementGaussPoints,
            self.maxNumberVolumeElementFamilies,
            self.numberVolumeElementFamilies,
            self.prestressFlag,
            self.numberVolumeElements,
            self.numberNodes)
            
        self.stateSize = self.elementSizeInfo[0]
        self.state0Size = self.elementSizeInfo[1]
        self.propertySize = self.elementSizeInfo[2]

        self.pointerToIen = None
	self.memorySize -= self.numberVolumeElementNodes* \
                           self.numberVolumeElements*self.intSize
        self.pointerToMat = None
	self.memorySize -= self.numberVolumeElements*self.intSize
        self.pointerToVolumeElementFamilyList = None
        self.memorySize -= 3*self.maxNumberVolumeElementFamilies*self.intSize
        self.pointerToIvftmp = None
        self.memorySize -= self.numberVolumeElementFamilies*self.intSize

        # Sort traction BC (unimplemented at present).


        # Sort split node entries.
        pylith3d.sort_split_nodes(
            self.pointerToNfault,
            self.pointerToIndxiel,
            self.numberSplitNodeEntries,
            self.numberVolumeElements)

        # Sort slippery node entries.
        pylith3d.sort_slip_nodes(
            self.pointerToNslip,
            self.pointerToIndxiel,
            self.numberSlipperyNodeEntries,
            self.numberVolumeElements)
            
	print ""
        print "Hello from pl3dsetup.sortmesh (end)!"

        return

        
    def sparsesetup(self, mesh):

        # This function sets up sparse matrix and associated storage.

        import pylith3d

	print ""
        print "Hello from pl3dsetup.sparsesetup (begin)!"
        print "Setting up sparse matrix storage:"
        
        # Initialize PETSc and set up logging
        import os
        wd = os.getcwd()
        pylith3d.PetscInitialize()
        os.chdir(wd) # Work-around MPI_Init() changing our directory!
        
        self.autoprestrStage, \
        self.elasticStage, \
        self.viscousStage, \
        self.iterateEvent = pylith3d.setupPETScLogging()

        # Initialize variables that are defined in this function.

        # Arrays to map element equation numbers to global
        self.pointerToLm = None
        self.pointerToLmx = None
        self.pointerToLmf = None

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

        # Temporary arrays that can be deallocated after use
        self.pointerToIndx = None
        self.pointerToLink = None
        self.pointerToNbrs = None

        # Localize global equation numbers in element index arrays.
        self.pointerToLm = pylith3d.allocateInt(
            self.numberDegreesFreedom*self.connectivitySize)
	self.memorySize += self.numberDegreesFreedom* \
	    self.connectivitySize*self.intSize
        self.pointerToLmx = pylith3d.allocateInt(
            self.numberDegreesFreedom*self.connectivitySize)
	self.memorySize += self.numberDegreesFreedom* \
	    self.connectivitySize*self.intSize
        self.pointerToLmf = pylith3d.allocateInt(
            self.connectivitySize)
	self.memorySize += self.connectivitySize*self.intSize

        pylith3d.local(
            self.pointerToId,
            self.numberNodes,
            self.pointerToIens,
            self.pointerToLm,
            self.numberVolumeElements,
            self.numberVolumeElementNodes)

        pylith3d.localf(
            self.pointerToIens,
            self.pointerToLmf,
            self.numberVolumeElements,
            self.pointerToNfault,
            self.numberSplitNodeEntries,
            self.numberVolumeElementNodes)

        pylith3d.localx(
            self.pointerToIdx,
            self.numberNodes,
            self.pointerToIens,
            self.pointerToLmx,
            self.numberVolumeElements,
            self.pointerToNslip,
            self.numberSlipperyNodeEntries,
            self.numberVolumeElementNodes)

        # Keeping this for now as it may be wanted for output
        # self.pointerToNslip = None
	# self.memorySize -= self.numberSlipDimensions* \
	#     self.numberSlipperyNodeEntries*self.intSize

        # Allocate and populate sparse matrix arrays.  Some of these are
        # temporary and are then deleted after use.
        self.workingArraySize = pylith3d.cmp_stiffsz(
            self.numberGlobalEquations,
            self.pointerToLm,
            self.pointerToLmx,
            self.numberVolumeElements,
            self.totalNumberSlipperyNodes,
            self.numberVolumeElementNodes)

        self.pointerToIndx = pylith3d.allocateInt(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.intSize
        self.pointerToLink = pylith3d.allocateInt(
            self.workingArraySize)
	self.memorySize += self.workingArraySize*self.intSize
        self.pointerToNbrs = pylith3d.allocateInt(
            self.workingArraySize)
	self.memorySize += self.workingArraySize*self.intSize

        self.stiffnessMatrixInfo = pylith3d.lnklst(
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

        self.A, self.rhs, self.sol = pylith3d.createPETScMat(mesh)
	self.memorySize += self.stiffnessMatrixSize*self.intSize

        self.stiffnessMatrixStats = pylith3d.makemsr(
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

	print ""
	print ""
        print "Sparse matrix information:"
	print ""
        print "numberGlobalEquations:     %i" % self.numberGlobalEquations
        print "workingArraySize:          %i" % self.workingArraySize
        print "stiffnessMatrixSize:       %i" % self.stiffnessTrueSize
        print "stiffnessOffDiagonalSize:  %i" % self.stiffnessOffDiagonalSize
        print "minimumNonzeroTermsPerRow: %i" % self.minimumNonzeroTermsPerRow
        print "maximumNonzeroTermsPerRow: %i" % self.maximumNonzeroTermsPerRow
        print "averageNonzeroTermsPerRow: %g" % self.averageNonzeroTermsPerRow
	print ""
        print "Hello from pl3dsetup.sparsesetup (end)!"

        return
        
    def allocateremaining(self):

        # This function allocates all remaining arrays that are needed for computations.
        
        import pylith3d

	print ""
        print "Hello from pl3dsetup.allocateremaining (begin)!"
        print "Allocating remaining storage:"
        
        # Initialize variables that are defined in this function.

        # Force/displacement vectors and list of force flags
        self.pointerToBextern = None
        self.pointerToBtraction = None
        self.pointerToBgravity = None
        self.pointerToBconcForce = None
        self.pointerToBintern = None
        self.pointerToBresid = None
        self.pointerToBwink = None
        self.pointerToBwinkx = None
        self.pointerToDispVec = None
        self.pointerToDprev = None
        self.listNforce = [0, 0, 0, 0, 0, 0, 0, 0]
        self.pointerToListArrayNforce = None

        # Body forces
        self.listGrav = [0.0, 0.0, 0.0]
        self.pointerToListArrayGrav = None

        # Displacement arrays and global dimensions
        self.pointerToD = None
        self.pointerToDeld = None
        self.pointerToDcur = None
        self.listNsysdat = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.pointerToListArrayNsysdat = None
        self.pointerToListArrayIddmat = None

        # Split node displacement arrays
        self.pointerToDfault = None
        self.pointerToTfault = None

        # Slippery node displacement arrays
        self.pointerToDx = None
        self.pointerToDeldx = None
        self.pointerToDxcur = None

        # Storage for element stiffness arrays
        self.pointerToS = None
        self.pointerToStemp = None

        # Element arrays and dimensions
        self.pointerToState = None
        self.pointerToDstate = None
        self.pointerToState0 = None
        self.pointerToDmat = None
        self.listNpar = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.pointerToListArrayNpar = None

        # Time step information
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

        # Tolerance information
        self.listRgiter = [0.0, 0.0, 0.0]
        self.pointerToListArrayRgiter = None
        

        # Create necessary lists and convert them to arrays
        self.listGrav = [
            self.gravityX.value,
            self.gravityY.value,
            self.gravityZ.value]
        self.pointerToListArrayGrav = pylith3d.doubleListToArray(
            self.listGrav)
	self.memorySize += 3*self.doubleSize
                             
        # Allocate memory for all additional arrays

        # Force vectors
        if self.numberTractionBc != 0:
            self.tractionFlag = 1
        if self.gravityX.value != 0.0 or self.gravityY.value != 0.0 or self.gravityZ.value != 0.0:
            self.gravityFlag = 1
        if self.numberConcForces != 0 or self.numberDifferentialForceEntries != 0:
            self.concForceFlag = 1
        if self.tractionFlag != 0 or self.gravityFlag != 0 or self.concForceFlag != 0:
            self.externFlag = 1
	if self.numberWinklerForces != 0:
	    self.winklerFlag = 1
	if self.numberSlipperyWinklerForces != 0:
	    self.slipperyWinklerFlag = 1

        self.pointerToBextern = pylith3d.allocateDouble(
            self.externFlag*self.numberGlobalEquations)
	self.memorySize += self.externFlag* \
                           self.numberGlobalEquations*self.doubleSize
        self.pointerToBtraction = pylith3d.allocateDouble(
            self.tractionFlag*self.numberGlobalEquations)
	self.memorySize += self.tractionFlag* \
                           self.numberGlobalEquations*self.doubleSize
        self.pointerToBgravity = pylith3d.allocateDouble(
            self.gravityFlag*self.numberGlobalEquations)
	self.memorySize += self.gravityFlag* \
                           self.numberGlobalEquations*self.doubleSize
        self.pointerToBconcForce = pylith3d.allocateDouble(
            self.concForceFlag*self.numberGlobalEquations)
	self.memorySize += self.concForceFlag* \
                           self.numberGlobalEquations*self.doubleSize
        self.pointerToBwink = pylith3d.allocateDouble(
            self.winklerFlag*self.numberGlobalEquations)
	self.memorySize += self.winklerFlag* \
                           self.numberGlobalEquations*self.doubleSize
        self.pointerToBwinkx = pylith3d.allocateDouble(
            self.slipperyWinklerFlag*self.numberGlobalEquations)
	self.memorySize += self.slipperyWinklerFlag* \
                           self.numberGlobalEquations*self.doubleSize
        self.pointerToBintern = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToBresid = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToDispVec = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToDprev = pylith3d.allocateDouble(
            self.usePreviousDisplacementFlag*self.numberGlobalEquations)
	self.memorySize += self.usePreviousDisplacementFlag* \
                           self.numberGlobalEquations*self.doubleSize
            
        # Displacement arrays
        self.pointerToD = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDeld = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDcur = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize

        # Slippery node arrays
        self.pointerToDx = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDeldx = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDxcur = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize

        # Split node arrays
        self.pointerToDfault = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSplitNodeEntries)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberSplitNodeEntries*self.doubleSize
        self.pointerToTfault = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSplitNodeEntries)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberSplitNodeEntries*self.doubleSize

        # Local stiffness matrix arrays
        self.pointerToS = pylith3d.allocateDouble(
            self.maxElementEquations*self.maxElementEquations)
	self.memorySize += self.maxElementEquations* \
	    self.maxElementEquations*self.doubleSize
        self.pointerToStemp = pylith3d.allocateDouble(
            self.maxElementEquations*self.maxElementEquations)
	self.memorySize += self.maxElementEquations* \
	    self.maxElementEquations*self.doubleSize

        # Element arrays
        self.pointerToState = pylith3d.allocateDouble(
            self.stateSize)
	self.memorySize += self.stateSize*self.doubleSize
        self.pointerToDstate = pylith3d.allocateDouble(
            self.stateSize)
	self.memorySize += self.stateSize*self.doubleSize
        self.pointerToDmat = pylith3d.allocateDouble(
            self.materialMatrixDimension*
            self.numberVolumeElementGaussPoints*
            self.numberVolumeElements)
	self.memorySize += self.materialMatrixDimension* \
	    self.numberVolumeElementGaussPoints* \
            self.numberVolumeElements*self.doubleSize
        self.pointerToListArrayIddmat = pylith3d.intListToArray( 
            self.listIddmat)
	self.memorySize += 36*self.intSize
        self.pointerToState0 = pylith3d.allocateDouble(
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
            self.winklerFlag,
            self.slipperyWinklerFlag,
            self.usePreviousDisplacementFlag]
        self.pointerToListArrayNforce = pylith3d.intListToArray(
            self.listNforce)
	self.memorySize += 8*self.intSize
           
        # ncodat array
        self.listNcodat = [
            self.analysisTypeInt,
            self.debuggingOutputInt]
        self.pointerToListArrayNcodat = pylith3d.intListToArray(
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
        self.pointerToListArrayNpar = pylith3d.intListToArray(
            self.listNpar)
	self.memorySize += 12*self.intSize

        # nprint array
        self.listNprint = [
            self.numberFullOutputs,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.ucdOutputInt]
        self.pointerToListArrayNprint = pylith3d.intListToArray(
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
        self.pointerToListArrayNsysdat = pylith3d.intListToArray(
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
        self.pointerToListArrayNunits = pylith3d.intListToArray(
            self.listNunits)
	self.memorySize += 6*self.intSize

        # nvisdat array
        self.listNvisdat = [
            self.numberCycles,
            self.numberTimeStepGroups,
            self.totalNumberTimeSteps,
            self.numberLoadHistories]
        self.pointerToListArrayNvisdat = pylith3d.intListToArray(
            self.listNvisdat)
	self.memorySize += 4*self.intSize
            
        # rgiter array
        self.listRgiter = [
            self.stressTolerance.value,
            self.minimumStrainPerturbation,
            self.initialStrainPerturbation]
        self.pointerToListArrayRgiter = pylith3d.doubleListToArray(
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
        self.pointerToListArrayRtimdat = pylith3d.doubleListToArray(
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
        self.pointerToListArrayNtimdat = pylith3d.intListToArray(
            self.listNtimdat)
        self.memorySize += 9*self.intSize

	print ""
        print "Hello from pl3dsetup.allocateremaining (end)!"

        return


    def meshwrite(self):

        # This function outputs mesh information.
        # In the near future, this needs to be broken into classes for
        # Ascii output, plot output, UCD output, etc.

        import pylith3d

	print ""
        print "Hello from pl3dsetup.meshwriteascii (begin)!"
        print "Outputting Ascii mesh information:"

        # Write out global parameters
        pylith3d.write_global_info(
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

        # Write out nodal coordinates
        pylith3d.write_coords(
            self.pointerToX,
            self.numberNodes,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.asciiOutputFile,
            self.plotOutputFile)

        # Write out nodal boundary condition info
        pylith3d.write_bc(
            self.pointerToBond,
            self.pointerToIbond,
            self.numberNodes,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        # Write out local coordinate rotations
        pylith3d.write_skew(
            self.pointerToSkew,
            self.numberRotationEntries,
            self.autoRotateSlipperyNodesInt,
            self.numberNodes,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        # Write stress computation and subiteration parameters.
        pylith3d.write_strscomp(
            self.stressTolerance.value,
            self.minimumStrainPerturbation,
            self.initialStrainPerturbation,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        pylith3d.write_subiter(
            self.usePreviousDisplacementFlag,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        # Write out time step information
        pylith3d.write_timdat(
            self.pointerToDelt,
            self.pointerToAlfa,
            self.pointerToUtol,
            self.pointerToFtol,
            self.pointerToEtol,
            self.pointerToTimes,
            self.pointerToMaxstp,
            self.pointerToMaxit,
            self.pointerToNtdinit,
            self.pointerToLgdef,
            self.pointerToItmax,
            self.numberTimeStepGroups,
            self.totalNumberTimeSteps,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        # Write out timesteps when full output is desired
        pylith3d.write_fuldat(
            self.pointerToIprint,
            self.numberFullOutputs,
            self.analysisTypeInt,
            self.numberCycles,
            self.totalNumberTimeSteps,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.asciiOutputFile,
            self.plotOutputFile)

        # Write out state variables desired for output
        pylith3d.write_stateout(
            self.pointerToIstatout,
            self.pointerToNstatout,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.asciiOutputFile,
            self.plotOutputFile)

        # Write out load history information and deallocate Times array
        pylith3d.write_hist(
            self.pointerToHistry,
            self.pointerToTimes,
            self.numberLoadHistories,
            self.totalNumberTimeSteps,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        self.Times = None
	self.memorySize -= (self.totalNumberTimeSteps+1)*self.doubleSize

        # Write element info
        pylith3d.write_element_info(
            self.numberVolumeElements,
            self.numberVolumeElementNodes,
	    self.numberVolumeElementGaussPoints,
            self.volumeElementType,
            self.quadratureOrderInt,
            self.prestressAutoComputeInt,
            self.prestressAutoChangeElasticPropsInt,
            self.prestressAutoComputePoisson,
            self.prestressAutoComputeYoungs.value,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        # Write element node array and deallocate Indxiel
        pylith3d.write_connect(
            self.pointerToIens,
            self.pointerToIvfamily,
            self.pointerToIndxiel,
            self.numberVolumeElementNodes,
	    self.numberVolumeElementGaussPoints,
            self.numberVolumeElements,
            self.volumeElementType,
            self.numberVolumeElementFamilies,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.asciiOutputFile,
            self.plotOutputFile)

        self.pointerToIndxiel = None
	self.memorySize -= self.numberVolumeElements*self.intSize

        # Write material properties
        pylith3d.write_props(
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

        # Write mesh info to UCD file, if requested
        if self.ucdOutputInt >= 0:
            pylith3d.write_ucd_mesh(
                self.pointerToX,
                self.numberNodes,
                self.pointerToIens,
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
                self.ucdOutputInt,
                self.ucdOutputRoot)

        # Write split node info
        pylith3d.write_split(
            self.pointerToFault,
            self.pointerToNfault,
            self.numberSplitNodeEntries,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.asciiOutputFile,
            self.plotOutputFile)

        # Write slippery node info
        pylith3d.write_slip(
            self.pointerToNslip,
            self.numberSlipperyNodeEntries,
            self.totalNumberSlipperyNodes,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.asciiOutputFile,
            self.plotOutputFile)

        # Write differential force info and deallocate Nslip
        pylith3d.write_diff(
            self.pointerToDiforc,
            self.pointerToNslip,
            self.pointerToIdhist,
            self.numberSlipperyNodeEntries,
            self.numberDifferentialForceEntries,
            self.numberNodes,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        self.pointerToNslip = None
        self.memorySize -= self.numberSlipDimensions* \
                           self.numberSlipperyNodeEntries*self.intSize

        # Write split nodes to plot file, if requested and deallocate Idftn
        pylith3d.write_split_plot(
            self.pointerToIdftn,
            self.totalNumberSplitNodes,
            self.f77PlotOutput,
            self.plotOutputInt,
            self.plotOutputFile)

        self.pointerToIdftn = None
	self.memorySize -= self.totalNumberSplitNodes*self.intSize

        # Write Winkler force info and deallocate definition arrays
        pylith3d.write_wink(
            self.pointerToWinkdef,
            self.pointerToIwinkdef,
            self.pointerToIwinkid,
            self.numberWinklerEntries,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        self.pointerToWinkdef = None
	self.memorySize -= self.numberDegreesFreedom* \
                           self.numberWinklerEntries*self.doubleSize
        self.pointerToIwinkdef = None
	self.memorySize -= self.numberDegreesFreedom* \
                           self.numberWinklerEntries*self.intSize

        # Write slippery node Winkler force info and deallocate definition arrays
        pylith3d.write_winkx(
            self.pointerToWinkxdef,
            self.pointerToIwinkxdef,
            self.pointerToIwinkxid,
            self.numberSlipperyWinklerEntries,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        self.pointerToWinkxdef = None
	self.memorySize -= self.numberDegreesFreedom* \
                           self.numberSlipperyWinklerEntries*self.doubleSize
        self.pointerToIwinkxdef = None
	self.memorySize -= self.numberDegreesFreedom* \
                           self.numberSlipperyWinklerEntries*self.intSize

        # Write sparse matrix info
        pylith3d.write_sparse_info(
            self.numberGlobalEquations,
            self.stiffnessMatrixSize,
            self.minimumNonzeroTermsPerRow,
            self.maximumNonzeroTermsPerRow,
            self.averageNonzeroTermsPerRow,
            self.asciiOutputInt,
            self.f77AsciiOutput,
            self.asciiOutputFile)

	print ""
        print "Hello from pl3dsetup.meshwrite (end)!"

        return

        
    def __init__(self):
        Component.__init__(self, "pl3dsetup", "setup")

	print ""
        print "Hello from pl3dsetup.__init__!"

        return



# version
# $Id: Pylith3d_setup.py,v 1.31 2005/06/24 20:22:03 willic3 Exp $

# End of file 
