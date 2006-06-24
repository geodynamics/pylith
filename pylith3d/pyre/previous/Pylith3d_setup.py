#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
#
#  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
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

        self.stressTolerance = pl3dscan.stressTolerance.value
        self.minimumStrainPerturbation = pl3dscan.minimumStrainPerturbation
        self.initialStrainPerturbation = pl3dscan.initialStrainPerturbation
        self.preconditionerType = pl3dscan.preconditionerType
        self.maxPcgIterations = pl3dscan.maxPcgIterations

        self.displacementAccuracyMult = pl3dscan.displacementAccuracyMult
        self.forceAccuracyMult = pl3dscan.forceAccuracyMult
        self.energyAccuracyMult = pl3dscan.energyAccuracyMult
        
        self.minDisplacementAccuracy = pl3dscan.minDisplacementAccuracy
        self.minForceAccuracy = pl3dscan.minForceAccuracy
        self.minEnergyAccuracy = pl3dscan.minEnergyAccuracy

        self.gravityX = pl3dscan.gravityX.value
        self.gravityY = pl3dscan.gravityY.value
        self.gravityZ = pl3dscan.gravityZ.value

        self.prestressAutoComputePoisson = pl3dscan.prestressAutoComputePoisson
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

                                                                   
        # Initialize and define some integer parameters based on string or logical parameters in python

        self.asciiOutputInt = 0
        if pl3dscan.inventory.asciiOutput == "none":
            self.asciiOutputInt = 0
        elif pl3dscan.inventory.asciiOutput == "echo":
            self.asciiOutputInt = 1
        else:
            self.asciiOutputInt = 2
            
        self.plotOutputInt = 0
        if pl3dscan.inventory.plotOutput == "ascii":
            self.plotOutputInt = 0
        else:
            self.plotOutputInt = 1
            
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

        self.preconditionerTypeInt = 0
        if pl3dscan.preconditionerType == "diagonalNoUpdate":
            self.preconditionerTypeInt = 1
        elif pl3dscan.preconditionerType == "gaussSeidelNoUpdate":
            self.preconditionerTypeInt = 2
        elif pl3dscan.preconditionerType == "diagonalUpdate":
            self.preconditionerTypeInt = 3
        elif pl3dscan.preconditionerType == "gaussSeidelUpdate":
            self.preconditionerTypeInt = 4
        else:
            self.preconditionerTypeInt = 1

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
        self.maxElementNodes = pl3dscan._maxElementNodes
        self.maxGaussPoints = pl3dscan._maxGaussPoints
        self.maxElementEquations = pl3dscan._maxElementEquations
        self.numberElementTypes = pl3dscan._numberElementTypes
        self.pointerToListArrayNumberElementNodesBase = pl3dscan._pointerToListArrayNumberElementNodesBase
        self.pointerToElementTypeInfo = pl3dscan._pointerToElementTypeInfo
        self.pointerToSh = pl3dscan._pointerToSh
        self.pointerToShj = pl3dscan._pointerToShj
        self.pointerToGauss = pl3dscan._pointerToGauss

        # Invariant parameters related to material model
        self.maxMaterialModels = pl3dscan._maxMaterialModels
        self.maxStateVariables = pl3dscan._maxStateVariables
        self.pointerToMaterialModelInfo = pl3dscan._pointerToMaterialModelInfo
        self.pointerToMaterialModelStates = pl3dscan._pointerToMaterialModelStates

        # Parameters derived from values in the inventory or the
        # category 2 parameters.
        self.analysisTypeInt = pl3dscan._analysisTypeInt
        self.quadratureOrderInt = pl3dscan._quadratureOrderInt
        self.prestressAutoComputeInt = pl3dscan._prestressAutoComputeInt


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
        self.propertyListSize = pl3dscan._propertyListSize
        self.pointerToListArrayPropertyList = pl3dscan._pointerToListArrayPropertyList
        self.pointerToMaterialInfo = pl3dscan._pointerToMaterialInfo

        self.numberElements = pl3dscan._numberElements
        self.connectivitySize = pl3dscan._connectivitySize

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

    def run(self):
        import pylith3d

        # This function performs a lot of memory allocation, and also bundles several
        # scalar values into lists.  It is assumed that the scalar values contained in
        # the lists will no longer be accessible (or useful) from within python once
        # they are passed as arrays.  The contents of the lists should be maintained as
        # global variables, however.

	print ""
        print "Hello from pl3dsetup.run (begin)!"
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
        self.pointerToListArrayWscal = pylith3d.doubleListToArray(
            self.listWscal)
	self.memorySize += 3*self.doubleSize

        self.listRmult = [
            self.displacementAccuracyMult,
            self.forceAccuracyMult,
            self.energyAccuracyMult]
        self.pointerToListArrayRmult = pylith3d.doubleListToArray(
            self.listRmult)
	self.memorySize += 3*self.doubleSize

        self.listRmin = [
            self.minDisplacementAccuracy,
            self.minForceAccuracy,
            self.minEnergyAccuracy]
        self.pointerToListArrayRmin = pylith3d.doubleListToArray(
            self.listRmin)
	self.memorySize += 3*self.doubleSize

        self.listGrav = [
            self.gravityX,
            self.gravityY,
            self.gravityZ]
        self.pointerToListArrayGrav = pylith3d.doubleListToArray(
            self.listGrav)
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

        # Node-based info (coordinates, ID arrays, displacement arrays, winkler arrays,
        #  BC, and skew BC).

        self.pointerToX = pylith3d.allocateDouble(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes* \
	    self.doubleSize
        self.pointerToId = pylith3d.allocateInt(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes* \
	    self.intSize
        self.pointerToIwink = pylith3d.allocateInt(
            2*self.numberWinklerForces)
	self.memorySize += 2*self.numberWinklerForces*self.intSize
        self.pointerToWink = pylith3d.allocateDouble(
            self.numberWinklerForces)
	self.memorySize += self.numberWinklerForces*self.doubleSize
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
        # For now, try allocating everything whether it is used or not.
        #        if self.numberRotationEntries != 0 or self.autoRotateSlipperyNodes:
        #            self.pointerToSkew = pylith3d.allocateDouble(
        #                self.numberSkewDimensions*self.numberNodes)
        self.pointerToSkew = pylith3d.allocateDouble(
            self.numberSkewDimensions*self.numberNodes)
	self.memorySize += self.numberSkewDimensions* \
	    self.numberNodes* \
	    self.doubleSize


        # try:

        pylith3d.read_coords(
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

        self.numberGlobalEquations = pylith3d.read_bc(
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

        pylith3d.read_wink(
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

        pylith3d.read_skew(
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
        pylith3d.write_strscomp(
            self.stressTolerance,
            self.minimumStrainPerturbation,
            self.initialStrainPerturbation,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)
        
        pylith3d.write_subiter(
            self.pointerToListArrayRmult,
            self.pointerToListArrayRmin,
            self.preconditionerTypeInt,
            self.maxPcgIterations,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)
                             
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
            4*self.maxStateVariables)
	self.memorySize += 4*self.maxStateVariables*self.intSize

        # try:
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
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.timeStepInputFile,
            self.asciiOutputFile)

        pylith3d.read_fuldat(
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

        pylith3d.read_stateout(
            self.pointerToIstatout,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.stateVariableInputFile,
            self.asciiOutputFile,
            self.plotOutputFile)

        pylith3d.read_hist(
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
        pylith3d.write_element_info(
            self.numberElements,
            self.quadratureOrderInt,
            self.prestressAutoComputeInt,
            self.prestressAutoComputePoisson,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        self.pointerToIen = pylith3d.allocateInt(
            self.connectivitySize)
	self.memorySize += self.connectivitySize*self.intSize
        self.pointerToInfiel = pylith3d.allocateInt(
            6*self.numberElements)
	self.memorySize += 6*self.numberElements*self.intSize
        self.pointerToIndmat = pylith3d.allocateInt(
            self.numberMaterials)
	self.memorySize += self.numberMaterials*self.intSize
        self.pointerToImgrp = pylith3d.allocateInt(
            self.numberMaterials)
	self.memorySize += self.numberMaterials*self.intSize
        self.pointerToMhist = pylith3d.allocateInt(
            self.propertyListSize)
	self.memorySize += self.propertyListSize*self.intSize

        # try:
        pylith3d.write_props(
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
        pylith3d.read_mathist(
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


        self.elementSizeInfo = pylith3d.read_connect(
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

        # print "Just after read_connect"
        # print "memorySize: %d" % self.memorySize
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
        self.pointerToNfault = pylith3d.allocateInt(
            3*self.numberSplitNodeEntries)
	self.memorySize += 3*self.numberSplitNodeEntries*self.intSize
        self.pointerToFault = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSplitNodeEntries)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberSplitNodeEntries*self.doubleSize
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
        # For now, try allocating everything whether it is used or not.
        #        if self.numberDifferentialForceEntries != 0:
        #            self.pointerToIdhist = pylith3d.allocateInt(
        #                self.numberNodes)
        #            self.pointerToDiforc = pylith3d.allocateDouble(
        #                self.numberDegreesFreedom*self.numberNodes)
        # try:
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
        # print "Just after read_traction"

        self.totalNumberSplitNodes = pylith3d.read_split(
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

        self.totalNumberSlipperyNodes = pylith3d.read_slip(
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

        pylith3d.read_diff(
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
        
        self.pointerToIdx = pylith3d.allocateInt(
            self.numberSpaceDimensions*self.numberNodes)
	self.memorySize += self.numberSpaceDimensions* \
	    self.numberNodes*self.intSize
        self.pointerToIdftn = pylith3d.allocateInt(
            self.totalNumberSplitNodes)
	self.memorySize += self.totalNumberSplitNodes*self.intSize
        self.pointerToWinkx = pylith3d.allocateDouble(
            self.numberSlipperyWinklerForces)
	self.memorySize += self.numberSlipperyWinklerForces*self.doubleSize
        self.pointerToIwinkx = pylith3d.allocateInt(
            2*self.numberSlipperyWinklerForces)
	self.memorySize += 2*self.numberSlipperyWinklerForces*self.intSize
        self.pointerToIdslp = pylith3d.allocateInt(
            self.numberNodes)
	self.memorySize += self.numberNodes*self.intSize
        self.pointerToIpslp = pylith3d.allocateInt(
            self.numberSlipNeighbors*self.totalNumberSlipperyNodes)
	self.memorySize += self.numberSlipNeighbors* \
	    self.totalNumberSlipperyNodes*self.intSize

        pylith3d.id_split(
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
                        
        self.numberGlobalEquations = pylith3d.adjid(
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
        pylith3d.read_winkx(
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
        
        self.pointerToLm = pylith3d.allocateInt(
            self.numberDegreesFreedom*self.connectivitySize)
	self.memorySize += self.numberDegreesFreedom* \
	    self.connectivitySize*self.intSize
        self.pointerToLmx = pylith3d.allocateInt(
            self.numberDegreesFreedom*self.connectivitySize)
	self.memorySize += self.numberDegreesFreedom* \
	    self.connectivitySize*self.intSize
        # For now, try allocating everything whether it is used or not.
        #        if self.numberSplitNodeEntries != 0:
        #            self.pointerToLmf = pylith3d.allocateInt(
        #                self.numberElementNodes*self.numberElements)
        self.pointerToLmf = pylith3d.allocateInt(
            self.connectivitySize)
	self.memorySize += self.connectivitySize*self.intSize

        pylith3d.local(
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

        pylith3d.localf(
            self.pointerToIen,
            self.pointerToLmf,
            self.pointerToInfiel,
            self.connectivitySize,
            self.numberElements,
            self.pointerToElementTypeInfo,
            self.pointerToNfault,
            self.numberSplitNodeEntries)
        # print "Just after localf"

        pylith3d.localx(
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
        self.workingArraySize = 5000*self.numberGlobalEquations       # This may need to be adjusted.
        self.pointerToIndx = pylith3d.allocateInt(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.intSize
        self.pointerToLink = pylith3d.allocateInt(
            self.workingArraySize)
	self.memorySize += self.workingArraySize*self.intSize
        self.pointerToNbrs = pylith3d.allocateInt(
            self.workingArraySize)
	self.memorySize += self.workingArraySize*self.intSize

        # print "Just before lnklst"
	# print "memorySize: %d" % self.memorySize
	# try:
        self.stiffnessMatrixInfo = pylith3d.lnklst(
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

        self.pointerToJa = pylith3d.allocateInt(
            self.stiffnessMatrixSize)
	self.memorySize += self.stiffnessMatrixSize*self.intSize

        # print "Just before makemsr"
	# print "stiffnessMatrixSize: %i" % self.stiffnessMatrixSize
	# print "stiffnessOffDiagonalSize: %i" % self.stiffnessOffDiagonalSize
	# print "memorySize: %d" % self.memorySize
        # pylith3d.makemsr(
            # self.pointerToJa,
            # self.pointerToIndx,
            # self.pointerToLink,
            # self.pointerToNbrs,
            # self.numberGlobalEquations,
            # self.stiffnessMatrixSize,
            # self.workingArraySize)

        self.stiffnessMatrixStats = pylith3d.makemsr(
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

        pylith3d.write_sparse_info(
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
        self.pointerToAlnz = pylith3d.allocateDouble(
            self.stiffnessMatrixSize)
	self.memorySize += self.stiffnessMatrixSize*self.doubleSize
        self.pointerToPcg = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToZcg = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize

        # Force vectors
        self.pointerToB = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToBtot = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToBres = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToGvec1 = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToGvec2 = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
        self.pointerToPvec = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize

        # Displacement arrays
        self.pointerToD = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDeld = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
	self.memorySize += self.numberDegreesFreedom* \
	    self.numberNodes*self.doubleSize
        self.pointerToDprev = pylith3d.allocateDouble(
            self.numberGlobalEquations)
	self.memorySize += self.numberGlobalEquations*self.doubleSize
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
            self.stateVariableDimension*self.stateSize)
	self.memorySize += self.stateVariableDimension* \
	    self.stateSize*self.doubleSize
        self.pointerToDstate = pylith3d.allocateDouble(
            self.stateVariableDimension*self.stateSize)
	self.memorySize += self.stateVariableDimension* \
	    self.stateSize*self.doubleSize
        self.pointerToDmat = pylith3d.allocateDouble(
            self.materialMatrixDimension*self.dmatSize)
	self.memorySize += self.materialMatrixDimension* \
	    self.dmatSize*self.doubleSize
        self.pointerToListArrayIddmat = pylith3d.intListToArray( 
            self.listIddmat)
	self.memorySize += 36*self.intSize
	# print "memorySize: %d" % self.memorySize


        # Create arrays from lists that will be needed for the solution

        # ncodat array
        self.listNcodat = [
            self.analysisTypeInt,
            self.debuggingOutputInt]
        self.pointerToListArrayNcodat = pylith3d.intListToArray(
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
        self.pointerToListArrayNpar = pylith3d.intListToArray(
            self.listNpar)
	self.memorySize += 11*self.intSize

        # nprint array
        self.listNprint = [
            self.numberFullOutputs,
            self.asciiOutputInt,
            self.plotOutputInt]
        self.pointerToListArrayNprint = pylith3d.intListToArray(
            self.listNprint)
	self.memorySize += 3*self.intSize

        # nsiter array
        self.listNsiter = [
            self.preconditionerTypeInt,
            self.maxPcgIterations]
        self.pointerToListArrayNsiter = pylith3d.intListToArray(
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
        self.pointerToListArrayNsysdat = pylith3d.intListToArray(
            self.listNsysdat)
	self.memorySize += 11*self.intSize

        # nunits array
        self.listNunits = [
            self.f77StandardInput,
            self.f77StandardOutput,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput]
        self.pointerToListArrayNunits = pylith3d.intListToArray(
            self.listNunits)
	self.memorySize += 5*self.intSize

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
            self.stressTolerance,
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
            self.prestressAutoComputePoisson]
        self.pointerToListArrayRtimdat = pylith3d.doubleListToArray(
            self.listRtimdat)
	self.memorySize += 3*self.doubleSize

	print ""
	print ""
        print "Sparse matrix information:"
        print "minimumNonzeroTermsPerRow: %i" % self.minimumNonzeroTermsPerRow
        print "maximumNonzeroTermsPerRow: %i" % self.maximumNonzeroTermsPerRow
        print "averageNonzeroTermsPerRow: %g" % self.averageNonzeroTermsPerRow

	# print "memorySize: %d" % self.memorySize
	print ""
	print ""
        print "Hello from pl3dsetup.run (end)!"
        return


    def __init__(self):
        Component.__init__(self, "pl3dsetup", "setup")

	print ""
        print "Hello from pl3dsetup.__init__!"

        return



# version
# $Id: Pylith3d_setup.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $

# End of file 
