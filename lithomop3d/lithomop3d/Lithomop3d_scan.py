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


class Lithomop3d_scan(Component):


    def __init__(self):
        Component.__init__(self, "lm3dscan", "scanner")

        from pyre.units.pressure import Pa
        from pyre.units.length import m
        from pyre.units.time import s

        print "Hello from lm3dscan.__init__ (begin)!"

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
#        self.prestressQuadrature = 1
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

        print "Hello from lm3dscan.__init__ (end)!"
        
        return

# derived or automatically-specified quantities (category 3)

    def _init(self, parent):

        from math import sqrt
        import pyre.units
        import lithomop3d
        import string
        import os

        from pyre.units.pressure import Pa
        from pyre.units.length import m
        from pyre.units.time import s

        uparser = pyre.units.parser()

        print "Hello from lm3dscan._init (begin)!"

        # Initial definitions of all parameters:
        self._geometryType = ""
        self._geometryTypeInt = 0
        self._numberSpaceDimensions = 0
        self._numberDegreesFreedom = 0
        self._numberStressComponents = 0
#        self._numberElementNodes = 0
        self._maxElementNodes = 0
        self._numberElementTypes = 0
##  This should probably be done purely in f77 (no lists)
##        self._numberElementNodes = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                    0, 0]
        self._materialMatrixSize = 0
        self._numberSkewDimensions = 0
        self._numberSlipDimensions = 0
        self._numberSlipNeighbors = 0
        self._numberTractionDirections = 0
#        self._numberMaterialProperties = 0
        self._numberMaterialModels = 0
        self._numberMaterialProperties = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self._numberStateVariables = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#        self._numberElementCoordinates = 0
##  This should probably be done purely in f77 (no lists)
##        self._numberElementCoordinates = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                          0, 0]
#        self._numberElementEquations = 0
##  This should probably be done purely in f77 (no lists)
##        self._numberElementEquations = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                        0, 0]
#        self._numberGaussPoints = 0
        self._maxGaussPoints = 0
##  This should probably be done purely in f77 (no lists)
##        self._numberGaussPoints = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
##                                   0, 0, 0, 0]
#        self._numberPrestressGaussPoints = 0
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
        self._numberMaterialTypes = 0
#        self._densityUnits = "densityUnitsInitial12345678901"
#        self._densityScaleFactor = 0.0
#        self._youngUnits = "youngUnitsInitial1234567890123"
#        self._youngScaleFactor = 0.0
#        self._viscosityCoefficientUnits = "viscosityCoefficientUnitsIniti"
#        self._viscosityCoefficientScaleFactor = 0.0
#        self._cohesionUnits = "cohesionUnitsInitial1234567890"
#        self._cohesionScaleFactor = 0.0
        self._materialInfo = [0, 0, 0, 0]
#        self._viscousFlagInt = 0
#        self._plasticFlagInt = 0
#        self._materialHistoryFlagInt = 0
        self._numberMaterialHistories = 0
        self._numberElements = 0
        self._connectivitySize = 0
        self._lmArraySize = 0
        self._materialPropertySize = 0
        self._stateVariableSize = 0
        self._numberPrestressEntries = 0
        self._numberTractionBc = 0
        self._tractionBcUnits = "tractionBcUnitsInitial12345678"
        self._tractionBcScaleFactor = 0.0
        self._numberSplitNodeEntries = 0
        self._numberSlipperyNodeEntries = 0
        self._numberDifferentialForceEntries = 0
        self._numberSlipperyWinklerEntries = 0
        self._numberSlipperyWinklerForces = 0
        self._analysisTypeInt = 0
        self._prestressAutoComputeInt = 0

        # First see if there is a keyword = value file, which may be used
        # to override parameters from __init__.
        
        if self.inventory.keywordEqualsValueFile == "None":
            self._keywordEqualsValueFile = self.inventory.fileRoot + ".keyval"
        else:
            self._keywordEqualsValueFile = self.inventory.keywordEqualsValueFile

        # print "global: ", globals()
        # print "local: ", locals()
        if os.path.isfile(self._keywordEqualsValueFile):
            execfile(self._keywordEqualsValueFile)

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
        # specified by using Lithomop3d.  The geometry type is only used for
        # f77 routines and not in pyre. An integer value is also defined
        # for use in f77 routines.
        # Define some integer values that are derived from string variables.

        self._geometryType = "3D"
        self._geometryTypeInt = 4

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

        # Parameters that vary depending on geometry type
        self._numberSpaceDimensions = 3
        self._numberDegreesFreedom = 3
        self._numberStressComponents = 6
        self._maxElementNodes = 20
        self._numberElementTypes = 62
        self._numberElementNodes = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0]
        self._materialMatrixSize = 0
        self._numberElementNodes = 8
        self._materialMatrixSize = 21
        self._numberSkewDimensions = 2
        self._numberSlipDimensions = 5
        self._numberSlipNeighbors = 4
        self._numberTractionDirections = 2
        self._numberMaterialProperties = 10
        self._listIddmat = [
            1, 2, 3, 4, 5, 6,
            2, 7, 8, 9,10,11,
            3, 8,12,13,14,15,
            4, 9,13,16,17,18,
            5,10,14,17,19,20,
            6,11,15,18,20,21]
        self._numberElementCoordinates = 24
        self._numberElementEquations = 24

# Easily-derivable parameters
        self._numberGaussPoints = quadratureOrder*quadratureOrder*quadratureOrder
        if prestressQuadrature == 1:
            self._numberPrestressGaussPoints = 1
        else:
            self._numberPrestressGaussPoints = self._numberGaussPoints

        if self._numberGaussPoints == 1:
            self._listGauss = [0.0, 0.0, 0.0, 8.0]
        elif self._numberGaussPoints == 8:
            gaussPoint = 1.0/self._doubleRoot3
            self._listGauss = [
                -gaussPoint, -gaussPoint,  gaussPoint, 1.0,
                 gaussPoint, -gaussPoint,  gaussPoint, 1.0,
                 gaussPoint,  gaussPoint,  gaussPoint, 1.0,
                -gaussPoint,  gaussPoint,  gaussPoint, 1.0,
                -gaussPoint, -gaussPoint, -gaussPoint, 1.0,
                 gaussPoint, -gaussPoint, -gaussPoint, 1.0,
                 gaussPoint,  gaussPoint, -gaussPoint, 1.0,
                -gaussPoint,  gaussPoint, -gaussPoint, 1.0]

# Parameters derived from the number of entries in a file.
        try:
            self._numberNodes = lithomop3d.scan_coords(
                self._numberSpaceDimensions,
                f77FileInput,
                self._coordinateUnits,
                self._coordinateInputFile)

            self._coordinateScaleString = \
                                        1.0*uparser.parse(string.strip(self._coordinateUnits))
            self._coordinateScaleFactor = \
                                        self._coordinateScaleString/pyre.units.SI.meter

            self._numberBcEntries = lithomop3d.scan_bc(
                self._numberDegreesFreedom,
                f77FileInput,
                self._displacementUnits,
                self._velocityUnits,
                self._forceUnits,
                self._bcInputFile)

            self._displacementScaleString = \
                                          1.0*uparser.parse(string.strip(self._displacementUnits))
            self._displacementScaleFactor = \
                                          self._displacementScaleString/pyre.units.SI.meter
            self._velocityScaleString = \
                                      1.0*uparser.parse(string.strip(self._velocityUnits))
            self._velocityScaleFactor = \
                                      self._velocityScaleString*pyre.units.SI.second/pyre.units.SI.meter
            self._forceScaleString = \
                                   1.0*uparser.parse(string.strip(self._forceUnits))
            self._forceScaleFactor = \
                                   self._forceScaleString/pyre.units.SI.newton

            self._numberWinklerEntries = lithomop3d.scan_wink(
                self._numberDegreesFreedom,
                self._numberWinklerForces,
                f77FileInput,
                self._winklerInputFile)

            self._numberRotationEntries = lithomop3d.scan_skew(
                self._numberSkewDimensions,
                f77FileInput,
                self._rotationUnits,
                self._rotationInputFile)

            if self._numberRotationEntries != 0:
                self._rotationScaleString = \
                                          1.0*uparser.parse(string.strip(self._rotationUnits))
                self._rotationScaleFactor = \
                                          self._rotationScaleString/pyre.units.SI.radian

            self._timeStepInfo = lithomop3d.scan_timdat(
                f77FileInput,
                self._timeUnits,
                self._timeStepInputFile)
            self._numberTimeStepGroups = self._timeStepInfo[0]
            self._totalNumberTimeSteps = self._timeStepInfo[1]

            self._timeScaleString = \
                                  1.0*uparser.parse(string.strip(self._timeUnits))
            self._timeScaleFactor = \
                                  self._timeScaleString/pyre.units.SI.second

            self._numberFullOutputs = lithomop3d.scan_fuldat(
                self._analysisTypeInt,
                self._totalNumberTimeSteps,
                f77FileInput,
                self._fullOutputInputFile)

            self._numberLoadHistories = lithomop3d.scan_hist(
                f77FileInput,
                self._loadHistoryInputFile)

            self._materialInfo = lithomop3d.scan_prop(
                self._numberMaterialProperties,
                f77FileInput,
                self._densityUnits,
                self._youngUnits,
                self._viscosityCoefficientUnits,
                self._cohesionUnits,
                self._materialPropertiesInputFile)

            self._numberMaterialTypes = self._materialInfo[0]
            self._viscousFlagInt = self._materialInfo[1]
            self._plasticFlagInt = self._materialInfo[2]
            self._materialHistoryFlagInt = self._materialInfo[3]

            self._densityScaleString = \
                                     1.0*uparser.parse(string.strip(self._densityUnits))
            self._densityScaleFactor = \
                                     self._densityScaleString*pyre.units.SI.meter**3/pyre.units.SI.kilogram

            self._youngScaleString = \
                                   1.0*uparser.parse(string.strip(self._youngUnits))
            self._youngScaleFactor = \
                                   self._youngScaleString/pyre.units.SI.pascal

            self._viscosityCoefficientScaleString = \
                                                  1.0*uparser.parse(string.strip\
                                                                    (self._viscosityCoefficientUnits))
            self._viscosityCoefficientScaleFactor = \
                                                  self._viscosityCoefficientScaleString/ \
                                                  (pyre.units.SI.pascal*pyre.units.SI.second)

            self._cohesionScaleString = \
                                      1.0*uparser.parse(string.strip(self._cohesionUnits))
            self._cohesionScaleFactor = \
                                      self._cohesionScaleString/pyre.units.SI.pascal

            self._numberElements = lithomop3d.scan_connect(
                self._numberElementNodes,
                f77FileInput,
                self._connectivityInputFile)

            self._numberPrestressEntries = lithomop3d.scan_prestr(
                self._numberStressComponents,
                self._numberPrestressGaussPoints,
                self._numberElements,
                self._prestressAutoComputeInt,
                f77FileInput,
                self._prestressInputFile)

            self._numberTractionBc = lithomop3d.scan_traction(
                self._numberElementNodes,
                self._numberTractionDirections,
                self._tractionBcUnits,
                f77FileInput,
                self._tractionInputFile)

            if self._numberTractionBc != 0:
                self._tractionBcScaleString = \
                                            1.0*uparser.parse(string.strip(self._tractionBcUnits))
                self._tractionBcScaleFactor = \
                                            self._tractionBcScaleString/pyre.units.SI.pascal

            self._numberSplitNodeEntries = lithomop3d.scan_split(
                self._numberDegreesFreedom,
                f77FileInput,
                self._splitNodeInputFile)

            self._numberSlipperyNodeEntries = lithomop3d.scan_slip(
                self._numberDegreesFreedom,
                f77FileInput,
                self._slipperyNodeInputFile)

            self._numberDifferentialForceEntries = lithomop3d.scan_diff(
                self._numberDegreesFreedom,
                self._numberSlipperyNodeEntries,
                f77FileInput,
                self._differentialForceInputFile)

            self._numberSlipperyWinklerEntries = lithomop3d.scan_winkx(
                self._numberDegreesFreedom,
                self._numberSlipperyWinklerForces,
                self._numberSlipperyNodeEntries,
                f77FileInput,
                self._slipperyWinklerInputFile)

        except IOError, error:
            print "Situation:", error
        except ValueError, error:
            print "Situation:", error
        except:
            print "Exception from Lithomop3d_scan!"
                
        print "Hello from lm3dscan._init (end)!"

        return


    class Inventory(Component.Inventory):

        import pyre.properties
#  I don't think this is being used in this section any more.
#        from pyre.units.pressure import pascal


        inventory = [
            pyre.properties.str(
                "title",
                default="LITHOMOP3D simulation"),

            pyre.properties.str(
                "fileRoot",
                default="../examples/lithomop3d_test1/lithomop3d_test1"),

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
                "loadHistoryInputFile",
                default="None"),

            pyre.properties.str(
                "materialPropertiesInputFile",
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
# specified by using Lithomop3d.
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
# $Id: Lithomop3d_scan.py,v 1.1 2004/04/14 21:22:47 willic3 Exp $

# End of file 
