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
        import pyre.units
        lm3dscan = scanner

        print "Hello from lm3dsetup.initialize (begin)!"

        # Parameters needed from Lithomop3d_scan.py.
        # These parameters, which may be modified in the run function, should all be
        # available to Lithomop3d_run.py.  Inventory items are not altered here and
        # may be accessed from Lithomop3d_scan.py.

        self.winklerScaleX = lm3dscan.winklerScaleX
        self.winklerScaleY = lm3dscan.winklerScaleY
        self.winklerScaleZ = lm3dscan.winklerScaleZ

        self.stressTolerance = lm3dscan.stressTolerance/pyre.units.SI.pascal
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

        self.quadratureOrder = lm3dscan.quadratureOrder

        self.gravityX = lm3dscan.gravityX*pyre.units.SI.second**2/pyre.units.SI.meter
        self.gravityY = lm3dscan.gravityY*pyre.units.SI.second**2/pyre.units.SI.meter
        self.gravityZ = lm3dscan.gravityZ*pyre.units.SI.second**2/pyre.units.SI.meter

        self.prestressAutoCompute = lm3dscan.prestressAutoCompute
        self.prestressQuadrature = lm3dscan.prestressQuadrature
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

                                                                   
        # Initialize and define some integer parameters based on string or logical parameters in python

        self.asciiOutput = lm3dscan.inventory.asciiOutput
        self.asciiOutputInt = 0
        if self.asciiOutput == "none":
            self.asciiOutputInt = 0
        elif self.asciiOutput == "echo":
            self.asciiOutputInt = 1
        else:
            self.asciiOutputInt = 2
            
        self.plotOutput = lm3dscan.inventory.plotOutput
        self.plotOutputInt = 0
        if self.plotOutput == "ascii":
            self.plotOutputInt = 0
        else:
            self.plotOutputInt = 1
            
        self.debuggingOutputInt = 0
        if lm3dscan.inventory.debuggingOutput:
            self.debuggingOutputInt = 1
        else:
            self.debuggingOutputInt = 0

        self.autoRotateSlipperyNodes = lm3dscan.inventory.autoRotateSlipperyNodes
        self.autoRotateSlipperyNodesInt = 0
        if self.autoRotateSlipperyNodes:
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
        
        # First get all filenames

        self.keywordEqualsValueFile = lm3dscan._keywordEqualsValueFile
        self.asciiOutputFile = lm3dscan._asciiOutputFile
        self.plotOutputFile = lm3dscan._plotOutputFile
        self.coordinateInputFile = lm3dscan._coordinateInputFile
        self.bcInputFile = lm3dscan._bcInputFile
        self.winklerInputFile = lm3dscan._winklerInputFile
        self.rotationInputFile = lm3dscan._rotationInputFile
        self.timeStepInputFile = lm3dscan._timeStepInputFile
        self.fullOutputInputFile = lm3dscan._fullOutputInputFile
        self.loadHistoryInputFile = lm3dscan._loadHistoryInputFile
        self.materialPropertiesInputFile = lm3dscan._materialPropertiesInputFile
        self.connectivityInputFile = lm3dscan._connectivityInputFile
        self.prestressInputFile = lm3dscan._prestressInputFile
        self.tractionInputFile = lm3dscan._tractionInputFile
        self.splitNodeInputFile = lm3dscan._splitNodeInputFile
        self.slipperyNodeInputFile = lm3dscan._slipperyNodeInputFile
        self.differentialForceInputFile = lm3dscan._differentialForceInputFile
        self.slipperyWinklerInputFile = lm3dscan._slipperyWinklerInputFile

        self.geometryType = lm3dscan._geometryType
        self.geometryTypeInt = lm3dscan._geometryTypeInt
        self.integerZero = lm3dscan._integerZero
        self.integerOne = lm3dscan._integerOne
        self.integerTwo = lm3dscan._integerTwo
        self.integerThree = lm3dscan._integerThree
        self.integerFour = lm3dscan._integerFour
        self.doubleZero = lm3dscan._doubleZero
        self.doubleOne = lm3dscan._doubleOne
        self.doubleTwo = lm3dscan._doubleTwo
        self.doubleThree = lm3dscan._doubleThree
        self.doubleFour = lm3dscan._doubleFour
        self.doubleThird = lm3dscan._doubleThird
        self.doubleRoot3 = lm3dscan._doubleRoot3

        self.numberSpaceDimensions = lm3dscan._numberSpaceDimensions
        self.numberDegreesFreedom = lm3dscan._numberDegreesFreedom
        self.numberStressComponents = lm3dscan._numberStressComponents
        self.numberElementNodes = lm3dscan._numberElementNodes
        self.materialMatrixSize = lm3dscan._materialMatrixSize
        self.numberSkewDimensions = lm3dscan._numberSkewDimensions
        self.numberSlipDimensions = lm3dscan._numberSlipDimensions
        self.numberSlipNeighbors = lm3dscan._numberSlipNeighbors
        self.numberTractionDirections = lm3dscan._numberTractionDirections
        self.numberMaterialProperties = lm3dscan._numberMaterialProperties
        self.listIddmat = lm3dscan._listIddmat
        self.numberElementCoordinates = lm3dscan._numberElementCoordinates
        self.numberElementEquations = lm3dscan._numberElementEquations
        self.numberGaussPoints = lm3dscan._numberGaussPoints
        self.numberPrestressGaussPoints = lm3dscan._numberPrestressGaussPoints
        self.listGauss = lm3dscan._listGauss

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

        self.numberMaterialTypes = lm3dscan._numberMaterialTypes
        self.densityScaleFactor = lm3dscan._densityScaleFactor
        self.youngScaleFactor = lm3dscan._youngScaleFactor
        self.viscosityCoefficientScaleFactor = lm3dscan._viscosityCoefficientScaleFactor
        self.cohesionScaleFactor = lm3dscan._cohesionScaleFactor
        self.viscousFlagInt = lm3dscan._viscousFlagInt
        self.plasticFlagInt = lm3dscan._plasticFlagInt
        self.materialHistoryFlagInt = lm3dscan._materialHistoryFlagInt

        self.numberElements = lm3dscan._numberElements
        self.numberPrestressEntries = lm3dscan._numberPrestressEntries

        self.numberTractionBc = lm3dscan._numberTractionBc
        self.tractionBcScaleFactor = lm3dscan._tractionBcScaleFactor

        self.numberSplitNodeEntries = lm3dscan._numberSplitNodeEntries
        self.numberSlipperyNodeEntries = lm3dscan._numberSlipperyNodeEntries
        self.numberDifferentialForceEntries = lm3dscan._numberDifferentialForceEntries
        self.numberSlipperyWinklerEntries = lm3dscan._numberSlipperyWinklerEntries
        self.numberSlipperyWinklerForces = lm3dscan._numberSlipperyWinklerForces

        self.analysisTypeInt = lm3dscan._analysisTypeInt
        self.prestressAutoComputeInt = lm3dscan._prestressAutoComputeInt

        print "Hello from lm3dsetup.initialize (end)!"

        return

    def run(self):
        import lithomop3d

        # This function performs a lot of memory allocation, and also bundles several
        # scalar values into lists.  It is assumed that the scalar values contained in
        # the lists will no longer be accessible (or useful) from within python once
        # they are passed as arrays.  The contents of the lists should be maintained as
        # global variables, however.

        print "Hello from lm3dsetup.run (begin)!"

        # Initialize parameters that are defined in this function

        self.numberGlobalEquations = 0
        self.currentNumberEquations = 0

        self.totalNumberSplitNodes = 0
        self.totalNumberSlipperyNodes = 0

        self.pointerToX = None
        self.pointerToId = None
        self.pointerToIdx = None
        self.pointerToD = None
        self.pointerToBond = None
        self.pointerToIbond = None
        self.pointerToWink = None
        self.pointerToIwink = None
        self.pointerToSkew = None

        self.pointerToMaxstp = None
        self.pointerToDelt = None
        self.pointerToAlfa = None
        self.pointerToMaxit = None
        self.pointerToMaxitc = None
        self.pointerToLgdef = None
        self.pointerToIbbar = None
        self.pointerToUtol = None
        self.pointerToFtol = None
        self.pointerToEtol = None
        self.pointerToItmax = None
        self.pointerToIprint = None
        self.pointerToTimes = None
        self.pointerToHistry = None

        self.pointerToProp = None
        self.pointerToIen = None
        self.pointerToMat = None
        self.pointerToInfin = None
        self.pointerToStn = None
        self.pointerToSt0 = None
        
        self.pointerToDcur = None
        self.pointerToDxcur = None
        self.pointerToEps = None
        self.pointerToDeps = None
        self.pointerToScur = None
        self.pointerToIelno = None
        self.pointerToIside = None
        self.pointerToIhistry = None
        self.pointerToPres = None
        self.pointerToPdir = None
        self.pointerToNfault = None
        self.pointerToFault = None
        self.pointerToDfault = None
        self.pointerToNslip = None
        self.pointerToIdhist = None
        self.pointerToDiforc = None

        self.pointerToIdftn = None
        self.pointerToWinkx = None
        self.pointerToIwinkx = None
        self.pointerToIdslp = None
        self.pointerToIpslp = None

        self.pointerToXtmp = None
        self.pointerToItmp = None
        self.pointerToItmp1 = None
        self.pointerToItmp2 = None

        self.pointerToLm = None
        self.pointerToLmx = None
        self.pointerToLmf = None

        self.numberMaterialDimensions = 0
        self.pointerToDmat = None

        self.pointerToDeld = None
        self.pointerToBeta = None
        self.pointerToDbeta = None
        self.pointerToBetb = None
        self.pointerToDbetb = None

        self.pointerToDx = None
        self.pointerToTfault = None
        self.pointerToDeldx = None

        self.workingArraySize = 0
        self.pointerToIndx = None
        self.pointerToLink = None
        self.pointerToNbrs = None

        self.stiffnessMatrixSize = 0
        self.stiffnessOffDiagonalSize = 0
        self.stiffnessMatrixInfo = [0, 0]
        self.pointerToJa = None
        self.minimumNonzeroTermsPerRow = 0
        self.maximumNonzeroTermsPerRow = 0
        self.averageNonzeroTermsPerRow = 0.0
        self.stiffnessMatrixStats = [0, 0, 0.0]
        self.pointerToAlnz = None

        self.pointerToB = None
        self.pointerToBtot = None
        self.pointerToBres = None
        self.pointerToGvec1 = None
        self.pointerToGvec2 = None
        self.pointerToPvec = None
        self.pointerToPcg = None
        self.pointerToZcg = None
        self.pointerToDprev = None

        self.pointerToS = None
        self.pointerToStemp = None


        # Make lists that are used as arrays in the f77 function calls below.

        self.listWscal = [
            self.winklerScaleX,
            self.winklerScaleY,
            self.winklerScaleZ]
        self.pointerToListArrayWscal = lithomop3d.doubleListToArray(
            self.listWscal)

        self.listRmult = [
            self.displacementAccuracyMult,
            self.forceAccuracyMult,
            self.energyAccuracyMult]
        self.pointerToListArrayRmult = lithomop3d.doubleListToArray(
            self.listRmult)

        self.listRmin = [
            self.minDisplacementAccuracy,
            self.minForceAccuracy,
            self.minEnergyAccuracy]
        self.pointerToListArrayRmin = lithomop3d.doubleListToArray(
            self.listRmin)

        self.listGrav = [
            self.gravityX,
            self.gravityY,
            self.gravityZ]
        self.pointerToListArrayGrav = lithomop3d.doubleListToArray(
            self.listGrav)

        self.listPrscal = [
            self.prestressScaleXx,
            self.prestressScaleYy,
            self.prestressScaleZz,
            self.prestressScaleXy,
            self.prestressScaleXz,
            self.prestressScaleYz]
        self.pointerToListArrayPrscal = lithomop3d.doubleListToArray(
            self.listPrscal)
                                  
        self.listWxscal = [
            self.winklerSlipScaleX,
            self.winklerSlipScaleY,
            self.winklerSlipScaleZ]
        self.pointerToListArrayWxscal = lithomop3d.doubleListToArray(
            self.listWxscal)

        # Write out global parameters
        
        lithomop3d.write_global_info(
            self.title,
            self.asciiOutputInt,
            self.plotOutputInt,
            self.geometryTypeInt,
            self.numberNodes,
            self.numberSpaceDimensions,
            self.numberDegreesFreedom,
            self.analysisTypeInt,
            self.debuggingOutputInt,
            self.numberStressComponents,
            self.numberElementNodes,
            self.f77AsciiOutput,
            self.f77PlotOutput,
            self.asciiOutputFile,
            self.plotOutputFile)

        # Coordinates, ID arrays, boundary conditions, winkler forces, and skew BC

        self.pointerToX = lithomop3d.allocateDouble(
            self.numberSpaceDimensions*self.numberNodes)
        self.pointerToId = lithomop3d.allocateInt(
            self.numberSpaceDimensions*self.numberNodes)
        self.pointerToIdx = lithomop3d.allocateInt(
            self.numberSpaceDimensions*self.numberNodes)
        self.pointerToD = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
        self.pointerToBond = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
        self.pointerToIbond = lithomop3d.allocateInt(
            self.numberDegreesFreedom*self.numberNodes)
        self.pointerToWink = lithomop3d.allocateDouble(
            self.numberWinklerForces)
        self.pointerToIwink = lithomop3d.allocateInt(
            2*self.numberWinklerForces)
        self.pointerToSkew = lithomop3d.allocateDouble(
            self.numberSkewDimensions*self.numberNodes)
        # For now, try allocating everything whether it is used or not.
        #        if self.numberRotationEntries != 0 or self.autoRotateSlipperyNodes:
        #            self.pointerToSkew = lithomop3d.allocateDouble(
        #                self.numberSkewDimensions*self.numberNodes)

        try:

            lithomop3d.read_coords(
                self.pointerToX,
                self.coordinateScaleFactor,
                self.numberSpaceDimensions,
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
                self.numberDegreesFreedom,
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
                self.numberDegreesFreedom,
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
                self.numberSkewDimensions,
                self.numberNodes,
                self.autoRotateSlipperyNodesInt,
                self.f77FileInput,
                self.f77AsciiOutput,
                self.asciiOutputInt,
                self.rotationInputFile,
                self.asciiOutputFile)

        except IOError, error:
            print "Situation:", error
        except ValueError, error:
            print "Situation:", error
        except:
            print "Exception in block between read_coords and read_skew!"

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
        self.pointerToMaxstp = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
        self.pointerToDelt = lithomop3d.allocateDouble(
            self.numberTimeStepGroups)
        self.pointerToAlfa = lithomop3d.allocateDouble(
            self.numberTimeStepGroups)
        self.pointerToMaxit = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
        self.pointerToMaxitc = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
        self.pointerToLgdef = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
        self.pointerToIbbar = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
        self.pointerToUtol = lithomop3d.allocateDouble(
            self.numberTimeStepGroups)
        self.pointerToFtol = lithomop3d.allocateDouble(
            self.numberTimeStepGroups)
        self.pointerToEtol = lithomop3d.allocateDouble(
            self.numberTimeStepGroups)
        self.pointerToItmax = lithomop3d.allocateInt(
            self.numberTimeStepGroups)
        self.pointerToIprint = lithomop3d.allocateInt(
            self.numberFullOutputs)
        self.pointerToTimes = lithomop3d.allocateDouble(
            self.totalNumberTimeSteps+1)
        self.pointerToHistry = lithomop3d.allocateDouble(
            (self.totalNumberTimeSteps+1)*self.numberLoadHistories)

        try:
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
                self.pointerToMaxitc,
                self.pointerToLgdef,
                self.pointerToIbbar,
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
            
        except IOError, error:
            print "Situation:", error
        except ValueError, error:
            print "Situation:", error
        except:
            print "Exception in block between read_timdat and read_hist!"

        # Allocate and read info on material properties, connectivities, and prestresses
        lithomop3d.write_element_info(
            self.numberElements,
            self.numberGaussPoints,
            self.prestressAutoComputeInt,
            self.prestressAutoComputePoisson,
            self.numberPrestressGaussPoints,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)

        self.pointerToProp = lithomop3d.allocateDouble(
            self.numberMaterialProperties*self.numberMaterialTypes)
        self.pointerToIen = lithomop3d.allocateInt(
            self.numberElementNodes*self.numberElements)
        self.pointerToMat = lithomop3d.allocateInt(
            self.numberElements)
        self.pointerToInfin = lithomop3d.allocateInt(
            self.numberElements)
        self.pointerToStn = lithomop3d.allocateDouble(
            self.numberStressComponents*self.numberGaussPoints*self.numberElements)
        self.pointerToSt0 = lithomop3d.allocateDouble(
            self.numberStressComponents*self.numberPrestressGaussPoints*self.numberElements)

        try:
            lithomop3d.read_prop(
                self.pointerToProp,
                self.pointerToListArrayGrav,
                self.densityScaleFactor,
                self.youngScaleFactor,
                self.viscosityCoefficientScaleFactor,
                self.cohesionScaleFactor,
                self.numberDegreesFreedom,
                self.numberMaterialProperties,
                self.numberMaterialTypes,
                self.asciiOutputInt,
                self.plotOutputInt,
                self.f77FileInput,
                self.f77AsciiOutput,
                self.f77PlotOutput,
                self.viscousFlagInt,
                self.plasticFlagInt,
                self.materialPropertiesInputFile,
                self.asciiOutputFile,
                self.plotOutputFile)

            lithomop3d.read_connect(
                self.pointerToIen,
                self.pointerToMat,
                self.pointerToInfin,
                self.numberElementNodes,
                self.numberElements,
                self.numberNodes,
                self.numberMaterialTypes,
                self.numberGaussPoints,
                self.f77FileInput,
                self.f77AsciiOutput,
                self.f77PlotOutput,
                self.asciiOutputInt,
                self.plotOutputInt,
                self.connectivityInputFile,
                self.asciiOutputFile,
                self.plotOutputFile)

            # print "Just after read_connect"
            lithomop3d.read_prestr(
                self.pointerToStn,
                self.pointerToSt0,
                self.pointerToListArrayPrscal,
                self.numberStressComponents,
                self.numberGaussPoints,
                self.numberPrestressGaussPoints,
                self.numberElements,
                self.numberPrestressEntries,
                self.prestressAutoComputeInt,
                self.asciiOutputInt,
                self.f77FileInput,
                self.f77AsciiOutput,
                self.prestressInputFile,
                self.asciiOutputFile)
            # print "Just after read_prestr"

        except IOError, error:
            print "Situation:", error
        except ValueError, error:
            print "Situation:", error
        except:
            print "Exception in block between read_prop and read_prestr!"

        # Allocate some displacement and stress/strain arrays, and read traction, split node,
        # and slippery node input files.
        self.pointerToDcur = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
        self.pointerToDxcur = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
        self.pointerToEps = lithomop3d.allocateDouble(
            self.numberStressComponents*self.numberGaussPoints*self.numberElements)
        self.pointerToDeps = lithomop3d.allocateDouble(
            self.numberStressComponents*self.numberGaussPoints*self.numberElements)
        self.pointerToScur = lithomop3d.allocateDouble(
            self.numberStressComponents*self.numberGaussPoints*self.numberElements)
        self.pointerToIelno = lithomop3d.allocateInt(
            self.numberTractionBc)
        self.pointerToIside = lithomop3d.allocateInt(
            self.numberTractionBc)
        self.pointerToIhistry = lithomop3d.allocateInt(
            self.numberTractionBc)
        self.pointerToPres = lithomop3d.allocateDouble(
            (self.numberElementNodes/2)*self.numberTractionBc)
        self.pointerToPdir = lithomop3d.allocateDouble(
            self.numberTractionDirections*self.numberTractionBc)
        self.pointerToNfault = lithomop3d.allocateInt(
            3*self.numberSplitNodeEntries)
        self.pointerToFault = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSplitNodeEntries)
        self.pointerToDfault = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSplitNodeEntries)
        self.pointerToNslip = lithomop3d.allocateInt(
            self.numberSlipDimensions*self.numberSlipperyNodeEntries)
        self.pointerToIdhist = lithomop3d.allocateInt(
            self.numberNodes)
        self.pointerToDiforc = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
        # For now, try allocating everything whether it is used or not.
        #        if self.numberDifferentialForceEntries != 0:
        #            self.pointerToIdhist = lithomop3d.allocateInt(
        #                self.numberNodes)
        #            self.pointerToDiforc = lithomop3d.allocateDouble(
        #                self.numberDegreesFreedom*self.numberNodes)
        try:
            lithomop3d.read_traction(
                self.pointerToPres,
                self.pointerToPdir,
                self.tractionBcScaleFactor,
                self.pointerToIelno,
                self.pointerToIside,
                self.pointerToIhistry,
                self.numberTractionBc,
                self.numberElementNodes,
                self.numberTractionDirections,
                self.f77FileInput,
                self.f77AsciiOutput,
                self.asciiOutputInt,
                self.tractionInputFile,
                self.asciiOutputFile)
            # print "Just after read_traction"

            self.totalNumberSplitNodes = lithomop3d.read_split(
                self.pointerToFault,
                self.pointerToNfault,
                self.numberSplitNodeEntries,
                self.numberDegreesFreedom,
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

            self.totalNumberSlipperyNodes = lithomop3d.read_slip(
                self.pointerToNslip,
                self.numberSlipperyNodeEntries,
                self.numberDegreesFreedom,
                self.numberNodes,
                self.autoRotateSlipperyNodesInt,
                self.numberSlipDimensions,
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
                self.numberDegreesFreedom,
                self.numberNodes,
                self.numberSlipDimensions,
                self.f77FileInput,
                self.f77AsciiOutput,
                self.asciiOutputInt,
                self.differentialForceInputFile,
                self.asciiOutputFile)
            # print "Just after read_diff"

        except IOError, error:
            print "Situation:", error
        except ValueError, error:
            print "Situation:", error
        except:
            print "Exception in block between read_traction and read_diff!"

        # Allocate id arrays for split and slippery nodes and differential winkler arrays, then
        # adjust equations for the presence of slippery nodes, and read in the differential
        # winkler forces.
        # Note:  all input is done after this section.  The remainder of this function allocates
        # memory for computations.
        
        self.pointerToIdftn = lithomop3d.allocateInt(
            self.totalNumberSplitNodes)
        self.pointerToWinkx = lithomop3d.allocateDouble(
            self.numberSlipperyWinklerForces)
        self.pointerToIwinkx = lithomop3d.allocateInt(
            2*self.numberSlipperyWinklerForces)
        self.pointerToIdslp = lithomop3d.allocateInt(
            self.numberNodes)
        self.pointerToIpslp = lithomop3d.allocateInt(
            self.numberSlipNeighbors*self.totalNumberSlipperyNodes)

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
                        
        self.numberGlobalEquations = lithomop3d.adjid(
            self.pointerToId,
            self.pointerToIdx,
            self.pointerToNslip,
            self.pointerToIdslp,
            self.numberSlipperyNodeEntries,
            self.numberDegreesFreedom,
            self.numberNodes,
            self.totalNumberSlipperyNodes,
            self.numberSlipDimensions,
            self.currentNumberEquations)
        # print "Just after adjid"

        if self.totalNumberSlipperyNodes != 0 and self.autoRotateSlipperyNodes:
            self.pointerToXtmp = lithomop3d.allocateDouble(
                self.totalNumberSlipperyNodes)
            self.pointerToItmp = lithomop3d.allocateInt(
                self.totalNumberSlipperyNodes)
            self.pointerToItmp1 = lithomop3d.allocateInt(
                self.totalNumberSlipperyNodes)
            self.pointerToItmp2 = lithomop3d.allocateInt(
                self.totalNumberSlipperyNodes)

            lithomop3d.nfind(
                self.pointerToX,
                self.pointerToXtmp,
                self.pointerToIdslp,
                self.pointerToIpslp,
                self.pointerToItmp,
                self.pointerToItmp1,
                self.pointerToItmp2,
                self.pointerToNslip,
                self.numberSpaceDimensions,
                self.numberDegreesFreedom,
                self.numberSlipDimensions,
                self.numberSlipNeighbors,
                self.numberSlipperyNodeEntries,
                self.totalNumberSlipperyNodes,
                self.numberNodes)
            # print "Just after nfind"

            self.pointerToXtmp = None
            self.pointerToItmp = None
            self.pointerToItmp1 = None
            self.pointerToItmp2 = None
                         
        try:
            lithomop3d.read_winkx(
                self.pointerToWinkx,
                self.pointerToListArrayWxscal,
                self.pointerToIwinkx,
                self.pointerToIdx,
                self.numberNodes,
                self.numberDegreesFreedom,
                self.numberSlipperyWinklerForces,
                self.numberSlipperyWinklerEntries,
                self.f77FileInput,
                self.f77AsciiOutput,
                self.asciiOutputInt,
                self.slipperyWinklerInputFile,
                self.asciiOutputFile)
            # print "Just after read_winkx"
                             
        except IOError, error:
            print "Situation:", error
        except ValueError, error:
            print "Situation:", error
        except:
            print "Exception in block read_winkx!"
        
        # Allocate slippery node array used for automatic rotation computation
        # and lm arrays, as well as material matrix array.  Then construct the
        # lm arrays.

        self.pointerToLm = lithomop3d.allocateInt(
            self.numberDegreesFreedom*self.numberElementNodes*self.numberElements)
        self.pointerToLmx = lithomop3d.allocateInt(
            self.numberDegreesFreedom*self.numberElementNodes*self.numberElements)

        self.pointerToLmf = lithomop3d.allocateInt(
            self.numberElementNodes*self.numberElements)
        # For now, try allocating everything whether it is used or not.
        #        if self.numberSplitNodeEntries != 0:
        #            self.pointerToLmf = lithomop3d.allocateInt(
        #                self.numberElementNodes*self.numberElements)
            
        self.numberMaterialDimensions = self.numberElements

        if self.viscousFlagInt == 0 and self.plasticFlagInt == 0:
            self.numberMaterialDimensions = self.numberMaterialTypes

        self.pointerToDmat = lithomop3d.allocateDouble(
            self.materialMatrixSize*self.numberGaussPoints*self.numberMaterialDimensions)

        lithomop3d.local(
            self.pointerToId,
            self.pointerToIen,
            self.pointerToLm,
            self.numberElementNodes,
            self.numberDegreesFreedom,
            self.numberElements,
            self.numberNodes)
        # print "Just after local"

        lithomop3d.localf(
            self.pointerToNfault,
            self.pointerToIen,
            self.pointerToLmf,
            self.numberSplitNodeEntries,
            self.numberElementNodes,
            self.numberElements)
        # print "Just after localf"

        lithomop3d.localx(
            self.pointerToIdx,
            self.pointerToIen,
            self.pointerToLmx,
            self.pointerToNslip,
            self.numberElementNodes,
            self.numberDegreesFreedom,
            self.numberSlipperyNodeEntries,
            self.numberElements,
            self.numberNodes,
            self.numberSlipDimensions)
        # print "Just after localx"

        # Allocate memory for viscous computation arrays.
        self.pointerToDeld = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)

        self.pointerToBeta = lithomop3d.allocateDouble(
            self.numberStressComponents*self.numberGaussPoints*self.numberElements)
        self.pointerToDbeta = lithomop3d.allocateDouble(
            self.numberStressComponents*self.numberGaussPoints*self.numberElements)
        # For now, try allocating everything whether it is used or not.
        #        if self.viscousFlagInt == 1:
        #            self.pointerToBeta = lithomop3d.allocateDouble(
        #                self.numberStressComponents*self.numberGaussPoints*self.numberElements)
        #            self.pointerToDbeta = lithomop3d.allocateDouble(
        #                self.numberStressComponents*self.numberGaussPoints*self.numberElements)

        self.pointerToBetb = lithomop3d.allocateDouble(
            self.numberStressComponents*self.numberGaussPoints*self.numberElements)
        self.pointerToDbetb = lithomop3d.allocateDouble(
            self.numberStressComponents*self.numberGaussPoints*self.numberElements)
        # For now, try allocating everything whether it is used or not.
        # if self.plasticFlagInt == 1:
        #     self.pointerToBetb = lithomop3d.allocateDouble(
        #         self.numberStressComponents*self.numberGaussPoints*self.numberElements)
        #     self.pointerToDbetb = lithomop3d.allocateDouble(
        #         self.numberStressComponents*self.numberGaussPoints*self.numberElements)

        self.pointerToDx = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)
        self.pointerToTfault = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberSplitNodeEntries)
        self.pointerToDeldx = lithomop3d.allocateDouble(
            self.numberDegreesFreedom*self.numberNodes)

        # Allocate and populate sparse matrix arrays.  Several of these arrays are temporary and may
        # be deallocated after use.
        self.workingArraySize = 100*self.numberGlobalEquations       # This may need to be adjusted.
        self.pointerToIndx = lithomop3d.allocateInt(
            self.numberGlobalEquations)
        self.pointerToLink = lithomop3d.allocateInt(
            self.workingArraySize)
        self.pointerToNbrs = lithomop3d.allocateInt(
            self.workingArraySize)

        # print "Just before lnklst"
        # print "self.numberGaussPoints: %i" % self.numberGaussPoints
        # print "self.pointerToLm: %i" % self.pointerToLm
        # print "self.pointerToLmx: %i" % self.pointerToLmx
        # print "self.pointerToIndx: %i" % self.pointerToIndx
        # print "self.pointerToLink: %i" % self.pointerToLink
        # print "self.pointerToNbrs: %i" % self.pointerToNbrs
        # print "self.numberElementEquations: %i" % self.numberElementEquations
        # print "self.numberGlobalEquations: %i" % self.numberGlobalEquations
        # print "self.numberElements: %i" % self.numberElements
        # print "self.workingArraySize: %i" % self.workingArraySize

        self.stiffnessMatrixInfo = lithomop3d.lnklst(
            self.pointerToLm,
            self.pointerToLmx,
            self.pointerToIndx,
            self.pointerToLink,
            self.pointerToNbrs,
            self.numberElementEquations,
            self.numberGlobalEquations,
            self.numberElements,
            self.workingArraySize)
        # print "Just after lnklst"
        self.stiffnessMatrixSize = self.stiffnessMatrixInfo[0]
        self.stiffnessOffDiagonalSize = self.stiffnessMatrixInfo[1]

        self.pointerToJa = lithomop3d.allocateInt(
            self.stiffnessMatrixSize)

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
        self.pointerToAlnz = lithomop3d.allocateDouble(
            self.stiffnessMatrixSize)

        # Allocate additional vectors needed for solution, and output sparse matrix info.
        self.pointerToB = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
        self.pointerToBtot = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
        self.pointerToBres = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
        self.pointerToGvec1 = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
        self.pointerToGvec2 = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
        self.pointerToPvec = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
        self.pointerToPcg = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
        self.pointerToZcg = lithomop3d.allocateDouble(
            self.numberGlobalEquations)

        self.pointerToDprev = lithomop3d.allocateDouble(
            self.numberGlobalEquations)
        # For now, allocate everythin whether it is used or not.
        # if self.preconditionerType == "diagonalUpdate" or self.preconditionerType == "gaussSeidelUpdate":
        # self.pointerToDprev = lithomop3d.allocateDouble(
        # self.numberGlobalEquations)

        self.pointerToS = lithomop3d.allocateDouble(
            self.numberElementEquations*self.numberElementEquations)
        self.pointerToStemp = lithomop3d.allocateDouble(
            self.numberElementEquations*self.numberElementEquations)

        # print "self.minimumNonzeroTermsPerRow: %i" % self.minimumNonzeroTermsPerRow
        # print "self.maximumNonzeroTermsPerRow: %i" % self.maximumNonzeroTermsPerRow
        # print "self.averageNonzeroTermsPerRow: %g" % self.averageNonzeroTermsPerRow

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
                       
        # Create arrays from lists that will be needed for the solution

        # gauss array
        self.pointerToListArrayGauss = lithomop3d.doubleListToArray(
            self.listGauss)

        # iddmat array
        self.pointerToListArrayIddmat = lithomop3d.intListToArray(
            self.listIddmat)

        # ncodat array
        self.listNcodat = [
            self.analysisTypeInt,
            self.debuggingOutputInt]
        self.pointerToListArrayNcodat = lithomop3d.intListToArray(
            self.listNcodat)
            
        # nconsts array
        self.listNconsts = [
            self.integerZero,
            self.integerOne,
            self.integerTwo,
            self.integerThree,
            self.integerFour]
        self.pointerToListArrayNconsts = lithomop3d.intListToArray(
            self.listNconsts)

        # ndimens array
        self.listNdimens = [
            self.numberSpaceDimensions,
            self.numberDegreesFreedom,
            self.numberStressComponents,
            self.numberElementNodes,
            self.geometryTypeInt,
            self.materialMatrixSize,
            self.numberSkewDimensions,
            self.numberSlipDimensions,
            self.numberSlipNeighbors,
            self.numberTractionDirections]
        self.pointerToListArrayNdimens = lithomop3d.intListToArray(
            self.listNdimens)

        # npar array
        self.listNpar = [
            self.numberElements,
            self.numberMaterialTypes,
            self.numberTractionBc,
            self.numberSlipperyNodeEntries,
            self.numberSplitNodeEntries,
            self.numberMaterialDimensions,
            self.prestressAutoComputeInt,
            self.numberPrestressGaussPoints,
            self.numberGaussPoints,
            self.numberDifferentialForceEntries]
        self.pointerToListArrayNpar = lithomop3d.intListToArray(
            self.listNpar)

        # nprint array
        self.listNprint = [
            self.numberFullOutputs,
            self.asciiOutputInt,
            self.plotOutputInt]
        self.pointerToListArrayNprint = lithomop3d.intListToArray(
            self.listNprint)

        # nsiter array
        self.listNsiter = [
            self.preconditionerTypeInt,
            self.maxPcgIterations]
        self.pointerToListArrayNsiter = lithomop3d.intListToArray(
            self.listNsiter)

        # nsysdat array
        self.listNsysdat = [
            self.numberNodes,
            self.numberGlobalEquations,
            self.stiffnessMatrixSize,
            self.numberElementEquations,
            self.numberElementCoordinates,
            self.numberRotationEntries,
            self.numberPrestressEntries,
            self.totalNumberSlipperyNodes,
            self.totalNumberSplitNodes,
            self.numberMaterialProperties,
            self.numberWinklerForces,
            self.numberSlipperyWinklerForces,
            self.autoRotateSlipperyNodesInt]
        self.pointerToListArrayNsysdat = lithomop3d.intListToArray(
            self.listNsysdat)

        # nunits array
        self.listNunits = [
            self.f77StandardInput,
            self.f77StandardOutput,
            self.f77FileInput,
            self.f77AsciiOutput,
            self.f77PlotOutput]
        self.pointerToListArrayNunits = lithomop3d.intListToArray(
            self.listNunits)

        # nvisdat array
        self.listNvisdat = [
            self.numberCycles,
            self.numberTimeStepGroups,
            self.totalNumberTimeSteps,
            self.numberLoadHistories]
        self.pointerToListArrayNvisdat = lithomop3d.intListToArray(
            self.listNvisdat)
            
        # rconsts array
        self.listRconsts = [
            self.doubleZero,
            self.doubleOne,
            self.doubleTwo,
            self.doubleThree,
            self.doubleFour,
            self.doubleThird,
            self.doubleRoot3]
        self.pointerToListArrayRconsts = lithomop3d.doubleListToArray(
            self.listRconsts)

        # rgiter array
        self.listRgiter = [
            self.stressTolerance,
            self.minimumStrainPerturbation,
            self.initialStrainPerturbation]
        self.pointerToListArrayRgiter = lithomop3d.doubleListToArray(
            self.listRgiter)

        # rmin array
        self.listRmin = [
            self.minDisplacementAccuracy,
            self.minForceAccuracy,
            self.minEnergyAccuracy]
        self.pointerToListArrayRmin = lithomop3d.doubleListToArray(
            self.listRmin)

        # rmult array
        self.listRmult = [
            self.displacementAccuracyMult,
            self.forceAccuracyMult,
            self.energyAccuracyMult]
        self.pointerToListArrayRmult = lithomop3d.doubleListToArray(
            self.listRmult)

        # rtimdat array
        self.currentTimeStepSize = 0.0
        self.currentAlfaParameter = 0.0
        self.listRtimdat = [
            self.currentTimeStepSize,
            self.currentAlfaParameter,
            self.prestressAutoComputePoisson]
        self.pointerToListArrayRtimdat = lithomop3d.doubleListToArray(
            self.listRtimdat)


        print "Hello from lm3dsetup.run (end)!"
        return


    def __init__(self):
        Component.__init__(self, "lm3dsetup", "setup")

        print "Hello from lm3dsetup.__init__!"

        return



# version
# $Id: Lithomop3d_setup.py,v 1.1 2004/04/14 21:22:47 willic3 Exp $

# End of file 
