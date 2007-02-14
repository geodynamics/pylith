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


from pyre.components.Component import Component
import os


class Scanner(Component):


    def __init__(self):
        Component.__init__(self, "pl3dscan", "scanner")

        import journal
        self.trace = journal.debug("pylith3d.trace")

        self.trace.log("Hello from pl3dscan.__init__!")
        
        self.macros = {}
        self.defineMacros({'rank': '0'})
        
        return


    def configureProperties(self, context):
        """set the values of all the properties in my inventory"""

        #
        # First, configure my properties as defined by the framework.
        #
        
        super(Scanner, self).configureProperties(context)
        
        #
        # Second, open input files.  Log I/O errors.
        #
        
        self._summaryIOError = self.CanNotOpenInputOutputFilesError()

        inputFile = self.inputFile
        inputFileStream = self.inputFileStream
        outputFile = self.outputFile
        macroString = self.macroString

        #                              open?   fatal?  label
        optional = self.IOFileCategory(True,   0,      "optional")
        unused   = self.IOFileCategory(False,  0,      "unused")
        required = self.IOFileCategory(True,   1,       None)
        
        Inventory = Scanner.Inventory

        analysisType = self.inventory.analysisType

        self._asciiOutputFile             = outputFile(Inventory.asciiOutputFile,            optional)
        self._plotOutputFile              = outputFile(Inventory.plotOutputFile,             optional)
        self._coordinateInputFile         = inputFile(Inventory.coordinateInputFile,         required)
        self._bcInputFile                 = inputFile(Inventory.bcInputFile,                 required)
        self._winklerInputFile            = inputFile(Inventory.winklerInputFile,            optional)
        self._rotationInputFile           = inputFile(Inventory.rotationInputFile,           optional)
        self._timeStepInputFile           = inputFile(Inventory.timeStepInputFile,           required)
        self._fullOutputInputFile         = inputFile(Inventory.fullOutputInputFile, analysisType == "fullSolution" and required or unused)
        self._stateVariableInputFile      = inputFile(Inventory.stateVariableInputFile,      required)
        self._loadHistoryInputFile        = inputFile(Inventory.loadHistoryInputFile,        optional)
        self._materialPropertiesInputFile = inputFile(Inventory.materialPropertiesInputFile, required)
        self._materialHistoryInputFile    = inputFile(Inventory.materialHistoryInputFile,    unused)
        self._connectivityInputFile       = inputFile(Inventory.connectivityInputFile,       required)
        self._prestressInputFile          = inputFile(Inventory.prestressInputFile,          unused)
        self._tractionInputFile           = inputFile(Inventory.tractionInputFile,           optional)
        self._splitNodeInputFile          = inputFile(Inventory.splitNodeInputFile,          optional)
        # Slippery nodes are not yet implemented in PyLith-0.8.
        self._slipperyNodeInputFile       = inputFile(Inventory.slipperyNodeInputFile,       unused)
        self._differentialForceInputFile  = inputFile(Inventory.differentialForceInputFile,  unused)
        self._slipperyWinklerInputFile    = inputFile(Inventory.slipperyWinklerInputFile,    unused)

        # The call to glob() is somewhat crude -- basically, determine
        # if any files might be in the way.
        self._ucdOutputRoot               = macroString(Inventory.ucdOutputRoot)
        from glob import glob
        ucdFiles = ([self._ucdOutputRoot + ".mesh.inp",
                     self._ucdOutputRoot + ".gmesh.inp",
                     self._ucdOutputRoot + ".mesh.time.prest.inp",
                     self._ucdOutputRoot + ".gmesh.time.prest.inp"]
                    + glob(self._ucdOutputRoot + ".mesh.time.[0-9][0-9][0-9][0-9][0-9].inp")
                    + glob(self._ucdOutputRoot + ".gmesh.time.[0-9][0-9][0-9][0-9][0-9].inp"))
        trait = Inventory.ucdOutputRoot
        for ucdFile in ucdFiles:
            try:
                stream = os.fdopen(os.open(ucdFile, os.O_WRONLY|os.O_CREAT|os.O_EXCL), "w")
            except (OSError, IOError), error:
                descriptor = self.inventory.getTraitDescriptor(trait.name)
                self._summaryIOError.openFailed(trait, descriptor,self._ucdOutputRoot + ".*mesh*.inp", error, required)
                break
            else:
                stream.close()
                os.remove(ucdFile)

        if self._summaryIOError.fatalIOErrors():
            self._summaryIOError.log(self._error)
            import sys
            sys.exit("%s: configuration error(s)" % self.name)

        #
        # Finally, decorate any missing traits with my name.
        #
        
        self._claim(context)

        return


    def _configure(self):
        super(Scanner, self)._configure()
        
        # get values for extra input (category 2)

        self.winklerScaleX = self.inventory.winklerScaleX
        self.winklerScaleY = self.inventory.winklerScaleY
        self.winklerScaleZ = self.inventory.winklerScaleZ
        
        self.stressTolerance = self.inventory.stressTolerance
        self.minimumStrainPerturbation = self.inventory.minimumStrainPerturbation
        self.initialStrainPerturbation = self.inventory.initialStrainPerturbation
        
        self.usePreviousDisplacementFlag = self.inventory.usePreviousDisplacementFlag
        
        self.quadratureOrder = self.inventory.quadratureOrder
        
        self.gravityX = self.inventory.gravityX
        self.gravityY = self.inventory.gravityY
        self.gravityZ = self.inventory.gravityZ
        
        self.prestressAutoCompute = self.inventory.prestressAutoCompute
        self.prestressAutoChangeElasticProps = self.inventory.prestressAutoChangeElasticProps
        self.prestressAutoComputePoisson = self.inventory.prestressAutoComputePoisson
        self.prestressAutoComputeYoungs = self.inventory.prestressAutoComputeYoungs
        
        self.prestressScaleXx = self.inventory.prestressScaleXx
        self.prestressScaleYy = self.inventory.prestressScaleYy
        self.prestressScaleZz = self.inventory.prestressScaleZz
        self.prestressScaleXy = self.inventory.prestressScaleXy
        self.prestressScaleXz = self.inventory.prestressScaleXz
        self.prestressScaleYz = self.inventory.prestressScaleYz
        
        self.winklerSlipScaleX = self.inventory.winklerSlipScaleX
        self.winklerSlipScaleY = self.inventory.winklerSlipScaleY
        self.winklerSlipScaleZ = self.inventory.winklerSlipScaleZ
        
        self.f77StandardInput = self.inventory.f77StandardInput
        self.f77StandardOutput = self.inventory.f77StandardOutput
        self.f77FileInput = self.inventory.f77FileInput
        self.f77AsciiOutput = self.inventory.f77AsciiOutput
        self.f77PlotOutput = self.inventory.f77PlotOutput
        self.f77UcdOutput = self.inventory.f77UcdOutput

        return


# derived or automatically-specified quantities (category 3)

    def initialize(self):

        from Materials import Materials
        import pyre.units
        import pylith3d
        import string
        from mpi import MPI_Comm_rank, MPI_COMM_WORLD
        
        rank = MPI_Comm_rank(MPI_COMM_WORLD)
        self.defineMacros({'rank': str(rank)})
        
        outputFile = self.outputFile
        Inventory = Scanner.Inventory
        optional = self.IOFileCategory(True,   0,      "optional")
        # Re-initialize these with the newly acquired rank.
        self._asciiOutputFile             = outputFile(Inventory.asciiOutputFile,            optional)
        self._plotOutputFile              = outputFile(Inventory.plotOutputFile,             optional)

        uparser = pyre.units.parser()
        matinfo = Materials()

        self.trace.log("Hello from pl3dscan.preinitialize (begin)!")
        
        self.trace.log("Scanning ascii files to determine dimensions:")

        # Define information needed from other functions:
        f77FileInput = self.f77FileInput
        prestressAutoCompute = self.prestressAutoCompute
        prestressAutoChangeElasticProps = self.prestressAutoChangeElasticProps
        quadratureOrder = self.quadratureOrder

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
        # Note:  eventually the variable below should disappear, and the
        # total size of the state variable array for each material model
        # should be used instead.  This means that all state variable
        # bookkeeping should be done within the material model routines.
        self._stateVariableDimension = 0
        self._materialMatrixDimension = 0
        self._numberSkewDimensions = 0
        self._numberSlipDimensions = 0
        self._numberSlipNeighbors = 0
        self._listIddmat = [0]

        # Invariant parameters related to element type
        self._numberElementTypes = 0
        self._numberElementTypesBase = 0
        self._numberElementNodesBase = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self._pointerToListArrayNumberElementNodesBase = None

        # Invariant parameters related to material model
        self._maxMaterialModels = 0
        self._maxStateVariables = 0
        self._maxState0Variables = 0
        self._pointerToMaterialModelInfo = None

        # Parameters derived from values in the inventory or the
        # category 2 parameters above.
        self._analysisTypeInt = 0
        self._pythonTimestep = 0
        self._prestressAutoComputeInt = 0
        self._prestressAutoChangeElasticPropsInt = 0
        self._pointerToSh = None
        self._pointerToShj = None
        self._pointerToGauss = None
        self._pointerToSh2d = None
        self._pointerToGauss2d = None

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

        self._volumeElementDimens = [0, 0, 0]
        self._numberVolumeElements = 0
        self._volumeElementType = 0
        self._numberVolumeElementFamilies = 0
        self._maxNumberVolumeElementFamilies = 0
        self._numberAllowedVolumeElementTypes = 0
        self._pointerToVolumeElementFamilyList = None

        self._numberPrestressEntries = 0

        self._numberTractionBc = 0
        self._tractionBcUnits = "tractionBcUnitsInitial12345678"
        self._tractionBcScaleFactor = 0.0
        self._tractionFlag = 0

        self._numberSplitNodeEntries = 0

        self._numberSlipperyNodeEntries = 0
        self._numberDifferentialForceEntries = 0
	self._slipperyWinklerInfo = [0, 0]
        self._numberSlipperyWinklerEntries = 0
        self._numberSlipperyWinklerForces = 0

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
        self._maxElementNodes2d = 4
        self._maxGaussPoints2d = 4
        self._numberElementTypes2d = 2
        self._numberElementTypesBase2d = 2
        self._numberElementNodesBase2d = [4, 3]
        self._pointerToListArrayNumberElementNodesBase2d = pylith3d.intListToArray(
            self._numberElementNodesBase2d)
	self._memorySize += self._numberElementTypesBase2d*self._intSize

        # Invariant parameters related to material model
        self._maxMaterialModels = 20
        self._maxStateVariables = 30
        self._maxState0Variables = 6
        self._pointerToMaterialModelInfo = pylith3d.allocateInt(
            6*self._maxMaterialModels)
	self._memorySize += 6*self._maxMaterialModels*self._intSize

        pylith3d.matmod_def(
            self._pointerToMaterialModelInfo)

        # Parameters derived from values in the inventory or the
        # category 2 parameters above.
        analysisType = self.inventory.analysisType
        analysisTypeMap = {
            "dataCheck":       0,
            "stiffnessFactor": 1,
            "elasticSolution": 2,
            "fullSolution":    3,
            }
        self._analysisTypeInt = analysisTypeMap[analysisType]

        if prestressAutoCompute:
            self._prestressAutoComputeInt = 1
        else:
            self._prestressAutoComputeInt = 0

        if prestressAutoChangeElasticProps:
            self._prestressAutoChangeElasticPropsInt = 1
        else:
            self._prestressAutoChangeElasticPropsInt = 0

        # Parameters derived from the number of entries in a file.
        self._numberNodes = pylith3d.scan_coords(
            f77FileInput,
            self._coordinateUnits,
            self._coordinateInputFile)

        self._coordinateScaleString = \
                                    uparser.parse(string.strip(self._coordinateUnits))
        self._coordinateScaleFactor = self._coordinateScaleString.value

        self._numberBcEntries = pylith3d.scan_bc(
            f77FileInput,
            self._displacementUnits,
            self._velocityUnits,
            self._forceUnits,
            self._bcInputFile)

        if self._numberBcEntries > 0:
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

        self._numberMaterials = matinfo.readprop(self._materialPropertiesInputFile)

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

        # At present, we assume that the number of element families is equal to
        # the number of material types used, since only one volume element type at a
        # time is allowed.
        self._numberAllowedVolumeElementTypes = 1
        self._maxNumberVolumeElementFamilies = self._numberAllowedVolumeElementTypes* \
                                               self._numberMaterials

        self._pointerToVolumeElementFamilyList = pylith3d.allocateInt(
            3*self._maxNumberVolumeElementFamilies)
        self._memorySize += 3*self._maxNumberVolumeElementFamilies*self._intSize

        self._volumeElementDimens = pylith3d.scan_connect(
            self._pointerToListArrayNumberElementNodesBase,
            self._pointerToMaterialModelInfo,
            self._pointerToListArrayMaterialModel,
            self._pointerToVolumeElementFamilyList,
            self._maxNumberVolumeElementFamilies,
	    self._numberMaterials,
            f77FileInput,
            self._connectivityInputFile)

        self._numberVolumeElements = self._volumeElementDimens[0]
        self._numberVolumeElementFamilies = self._volumeElementDimens[1]
        self._volumeElementType = self._volumeElementDimens[2]

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

        self._numberTractionBc = pylith3d.scan_tractions(
            self._maxElementNodes2d,
            f77FileInput,
            self._tractionBcUnits,
            self._tractionInputFile)

        if self._numberTractionBc != 0:
            self._tractionBcScaleString = \
                                        uparser.parse(string.strip(self._tractionBcUnits))
            self._tractionBcScaleFactor = self._tractionBcScaleString.value
            self._tractionFlag = 1

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

        self.trace.log("Hello from pl3dscan.preinitialize (end)!")

        return


    class CanNotOpenInputOutputFilesError(Exception):
        
        def __init__(self):
            self._ioErrors = {}
            self._ioErrorProtos = {}
            self._fatalIOErrors = 0
        
        def openFailed(self, trait, descriptor, value, error, category):
            """Open failed for an I/O file property."""
            errno = error[0]
            if not self._ioErrors.has_key(errno):
                self._ioErrors[errno] = {}
                proto = IOError(error[0], error[1]) # omit filename
                self._ioErrorProtos[errno] = proto
            from copy import copy
            descriptor = copy(descriptor)
            descriptor.origValue = descriptor.value
            descriptor.value = value
            self._ioErrors[errno][trait.name] = (error, descriptor, category)
            self._fatalIOErrors = self._fatalIOErrors + category.fatalPoints
            return

        def fatalIOErrors(self): return self._fatalIOErrors

        def log(self, diag):
            from StringIO import StringIO
            diag.line("error(s) opening input/output file(s)")
            diag.line("")
            errnos = self._ioErrors.keys()
            errnos.sort()
            cw = [4, 4, 30, 15, 10, 10] # column widths
            ch = ("", "", "property", "value", "from", "") # column headers
            for errno in errnos:
                propertyNames = self._ioErrors[errno].keys()
                for name in propertyNames:
                    error, descriptor, category = self._ioErrors[errno][name]
                    valueLen = len(descriptor.value)
                    if valueLen > cw[3]:
                        cw[3] = valueLen
            for errno in errnos:
                stream = StringIO()
                print >> stream, "".ljust(cw[0]), self._ioErrorProtos[errno],
                diag.line(stream.getvalue())
                stream = StringIO()
                for column in xrange(0, len(ch)):
                    print >> stream, ch[column].ljust(cw[column]),
                diag.line(stream.getvalue())
                stream = StringIO()
                for column in xrange(0, len(ch)):
                    print >> stream, ("-" * len(ch[column])).ljust(cw[column]),
                diag.line(stream.getvalue())
                propertyNames = self._ioErrors[errno].keys()
                propertyNames.sort()
                for name in propertyNames:
                    error, descriptor, category = self._ioErrors[errno][name]
                    stream = StringIO()
                    print >> stream, \
                          "".ljust(cw[0]), \
                          "".ljust(cw[1]), \
                          name.ljust(cw[2]), \
                          descriptor.value.ljust(cw[3]), \
                          str(descriptor.locator).ljust(cw[4]), \
                          (category.label and ("(%s)" % category.label).ljust(cw[5]) or ""),
                    diag.line(stream.getvalue())
                diag.line("")
            diag.log()
            return

    class IOFileCategory(object):
        def __init__(self, tryOpen, fatalPoints, label):
            self.tryOpen = tryOpen
            self.fatalPoints = fatalPoints
            self.label = label

    def defineMacros(self, macros):
        self.macros.update(macros)
    
    def macroString(self, trait):
        from pyre.util import expandMacros
        descriptor = self.inventory.getTraitDescriptor(trait.name)
        return expandMacros(descriptor.value, self.macros)

    def ioFileStream(self, trait, flags, mode, category, opener=None):
        value = self.macroString(trait)
        stream = None
        if category.tryOpen:
            try:
                if opener is None:
                    stream = os.fdopen(os.open(value, flags), mode)
                else:
                    stream = opener(value, flags, mode)
            except (OSError, IOError), error:
                descriptor = self.inventory.getTraitDescriptor(trait.name)
                self._summaryIOError.openFailed(trait, descriptor, value, error, category)
        return value, stream

    def inputFile(self, trait, category):
        value, stream = self.ioFileStream(trait,os. O_RDONLY, "r", category)
        if stream is not None:
            stream.close()
        return value
    
    def inputFileStream(self, trait, category, opener=None):
        return self.ioFileStream(trait, os.O_RDONLY, "r", category, opener=opener)[1]
    
    def outputFile(self, trait, category):
        value, stream = self.ioFileStream(trait, os.O_WRONLY|os.O_CREAT|os.O_EXCL, "w", category)
        if stream is not None:
            stream.close()
            os.remove(value)
        return value

    class Inventory(Component.Inventory):

        import pyre.inventory
        MacroString = pyre.inventory.str
        OutputFile = pyre.inventory.str
        InputFile = pyre.inventory.str

        # Title
        title = pyre.inventory.str("title",
                                   default="PyLith-0.8 Simulation")
        title.meta['tip'] = "Title for this simulation"

        # Output filenames (all are optional).
        asciiOutputFile = OutputFile("asciiOutputFile",default="${outputFileRoot}.ascii")
        asciiOutputFile.meta['tip'] = "Pathname for ascii output file (overrides default from outputFileRoot)."

        plotOutputFile = OutputFile("plotOutputFile",default="${outputFileRoot}.plot")
        plotOutputFile.meta['tip'] = "Pathname for plot output file (overrides default from outputFileRoot)."

        ucdOutputRoot = MacroString("ucdOutputRoot",default="${outputFileRoot}")
        ucdOutputRoot.meta['tip'] = "Base name for UCD output files (overrides default from outputFileRoot)."

        # Required input files.
        coordinateInputFile = InputFile("coordinateInputFile",default="${inputFileRoot}.coord")
        coordinateInputFile.meta['tip'] = "Pathname for coordinate input file (overrides default from inputFileRoot)."

        bcInputFile = InputFile("bcInputFile",default="${inputFileRoot}.bc")
        bcInputFile.meta['tip'] = "Pathname for boundary condition input file (overrides default from inputFileRoot)."

        timeStepInputFile = InputFile("timeStepInputFile",default="${inputFileRoot}.time")
        timeStepInputFile.meta['tip'] = "Pathname for time step definitions input file (overrides default from inputFileRoot)."

        stateVariableInputFile = InputFile("stateVariableInputFile",default="${inputFileRoot}.statevar")
        stateVariableInputFile.meta['tip'] = "Pathname for file defining which state variables to output (overrides default from inputFileRoot)."

        materialPropertiesInputFile = InputFile("materialPropertiesInputFile",default="${inputFileRoot}.prop")
        materialPropertiesInputFile.meta['tip'] = "Pathname for file defining material properties (overrides default from inputFileRoot)."

        connectivityInputFile = InputFile("connectivityInputFile",default="${inputFileRoot}.connect")
        connectivityInputFile.meta['tip'] = "Pathname for connectivity input file (overrides default from inputFileRoot)."

        # This file is only required for time-dependent problems.
        fullOutputInputFile = InputFile("fullOutputInputFile",default="${inputFileRoot}.fuldat")
        fullOutputInputFile.meta['tip'] = "Pathname for file defining when to provide output (overrides default from inputFileRoot)."

        # Optional input files.
        winklerInputFile = InputFile("winklerInputFile",default="${inputFileRoot}.wink")
        winklerInputFile.meta['tip'] = "Pathname for Winkler force input file (overrides default from inputFileRoot)."

        rotationInputFile = InputFile("rotationInputFile",default="${inputFileRoot}.skew")
        rotationInputFile.meta['tip'] = "Pathname for skew rotations input file (overrides default from inputFileRoot)."

        loadHistoryInputFile = InputFile("loadHistoryInputFile",default="${inputFileRoot}.hist")
        loadHistoryInputFile.meta['tip'] = "Pathname for file defining load histories (overrides default from inputFileRoot)."

        splitNodeInputFile = InputFile("splitNodeInputFile",default="${inputFileRoot}.split")
        splitNodeInputFile.meta['tip'] = "Pathname for split node input file (overrides default from inputFileRoot)."

        # Unused input files.
        materialHistoryInputFile = InputFile("materialHistoryInputFile",default="${inputFileRoot}.mhist")
        materialHistoryInputFile.meta['tip'] = "Pathname for file defining material histories (overrides default from inputFileRoot -- presently unused)."

        prestressInputFile = InputFile("prestressInputFile",default="${inputFileRoot}.prestr")
        prestressInputFile.meta['tip'] = "Pathname for prestress input file (overrides default from inputFileRoot -- presently unused)."

        tractionInputFile = InputFile("tractionInputFile",default="${inputFileRoot}.traction")
        tractionInputFile.meta['tip'] = "Pathname for traction BC input file (overrides default from inputFileRoot)."

        slipperyNodeInputFile = InputFile("slipperyNodeInputFile",default="${inputFileRoot}.slip")
        slipperyNodeInputFile.meta['tip'] = "Pathname for slippery node input file (overrides default from inputFileRoot -- presently unused)."

        differentialForceInputFile = InputFile("differentialForceInputFile",default="${inputFileRoot}.diff")
        differentialForceInputFile.meta['tip'] = "Pathname for file defining slippery node differential forces (overrides default from inputFileRoot -- presently unused)."

        slipperyWinklerInputFile = InputFile("slipperyWinklerInputFile",default="${inputFileRoot}.winkx")
        slipperyWinklerInputFile.meta['tip'] = "Pathname for file defining slippery node Winkler forces (overrides default from inputFileRoot -- presently unused)."

        # Output option flags.
        asciiOutput = pyre.inventory.str("asciiOutput",default="echo")
        asciiOutput.validator = pyre.inventory.choice(["none","echo","full"])
        asciiOutput.meta['tip'] = "Type of ascii output desired (none, echo, full)."

        plotOutput = pyre.inventory.str("plotOutput",default="none")
        plotOutput.validator = pyre.inventory.choice(["none","ascii","binary"])
        plotOutput.meta['tip'] = "Type of plot output desired (none, ascii, binary)."

        ucdOutput = pyre.inventory.str("ucdOutput",default=None)
        ucdOutput.validator = pyre.inventory.choice(["none","ascii","binary"])
        ucdOutput.meta['tip'] = "Type of UCD output desired (none, ascii, binary)."

        # Additional option flags.
        analysisType = pyre.inventory.str("analysisType",default="fullSolution")
        analysisType.validator = pyre.inventory.choice(["dataCheck","stiffnessFactor",
                                                        "elasticSolution","fullSolution"])
        analysisType.meta['tip'] = "Type of analysis (dataCheck, stiffnessFactor, elasticSolution, fullSolution)."

        pythonTimestep = pyre.inventory.bool("pythonTimestep",default=False)
        pythonTimestep.meta['tip'] = "Whether to use python timestepping loop (enables VTK output for time-dependent solution)."

        debuggingOutput = pyre.inventory.bool("debuggingOutput",default=False)
        debuggingOutput.meta['tip'] = "Whether to produce debugging output."

        numberCycles = pyre.inventory.int("numberCycles",default=1)
        numberCycles.meta['tip'] = "Number of cycles of the given timestep definitions to perform (default=1)."

        interpolateMesh = pyre.inventory.bool("interpolateMesh",default=False)
        interpolateMesh.meta['tip'] = "Create intermediate mesh entities, such as edges and faces."

        partitioner = pyre.inventory.str("partitioner",default="chaco")
        partitioner.validator = pyre.inventory.choice(["chaco","parmetis"])
        partitioner.meta['tip'] = "Partitioner (chaco, parmetis)."

        # Unused option flags.
        autoRotateSlipperyNodes = pyre.inventory.bool("autoRotateSlipperyNodes",default=True)
        autoRotateSlipperyNodes.meta['tip'] = "Whether to performa automatic rotation for slippery nodes (presently unused)."

        #
        # category 2 parameters traditionally placed in *.keyval files
        #

        from pyre.units.pressure import Pa
        from pyre.units.length import m
        from pyre.units.time import s

        winklerScaleX = pyre.inventory.float("winklerScaleX", default=1.0)
        winklerScaleY = pyre.inventory.float("winklerScaleY", default=1.0)
        winklerScaleZ = pyre.inventory.float("winklerScaleZ", default=1.0)

        stressTolerance = pyre.inventory.dimensional("stressTolerance", default=1.0e-12*Pa)
        minimumStrainPerturbation = pyre.inventory.float("minimumStrainPerturbation", default=1.0e-7)
        initialStrainPerturbation = pyre.inventory.float("initialStrainPerturbation", default=1.0e-1)

        usePreviousDisplacementFlag = pyre.inventory.int("usePreviousDisplacementFlag", default=0)

        quadratureOrder = pyre.inventory.str("quadratureOrder", default="Full")
        quadratureOrder.validator = pyre.inventory.choice(["Full", "Reduced", "Selective"])

        gravityX = pyre.inventory.dimensional("gravityX", default=0.0*m/(s*s))
        gravityY = pyre.inventory.dimensional("gravityY", default=0.0*m/(s*s))
        gravityZ = pyre.inventory.dimensional("gravityZ", default=0.0*m/(s*s))

        prestressAutoCompute = pyre.inventory.bool("prestressAutoCompute", default=False)
        prestressAutoChangeElasticProps = pyre.inventory.bool("prestressAutoChangeElasticProps", default=False)
        prestressAutoComputePoisson = pyre.inventory.float("prestressAutoComputePoisson", default=0.49)
        prestressAutoComputeYoungs = pyre.inventory.dimensional("prestressAutoComputeYoungs", default=1.0e30*Pa)

        prestressScaleXx = pyre.inventory.float("prestressScaleXx", default=1.0)
        prestressScaleYy = pyre.inventory.float("prestressScaleYy", default=1.0)
        prestressScaleZz = pyre.inventory.float("prestressScaleZz", default=1.0)
        prestressScaleXy = pyre.inventory.float("prestressScaleXy", default=1.0)
        prestressScaleXz = pyre.inventory.float("prestressScaleXz", default=1.0)
        prestressScaleYz = pyre.inventory.float("prestressScaleYz", default=1.0)

        winklerSlipScaleX = pyre.inventory.float("winklerSlipScaleX", default=1.0)
        winklerSlipScaleY = pyre.inventory.float("winklerSlipScaleY", default=1.0)
        winklerSlipScaleZ = pyre.inventory.float("winklerSlipScaleZ", default=1.0)

        f77StandardInput = pyre.inventory.int("f77StandardInput", default=5)
        f77StandardOutput = pyre.inventory.int("f77StandardOutput", default=6)
        f77FileInput = pyre.inventory.int("f77FileInput", default=10)
        f77AsciiOutput = pyre.inventory.int("f77AsciiOutput", default=11)
        f77PlotOutput = pyre.inventory.int("f77PlotOutput", default=12)
        f77UcdOutput = pyre.inventory.int("f77UcdOutput", default=13)



from pyre.applications import ComponentHarnessAdapter
#from pyre.applications.ComponentHarness import ComponentHarness as ComponentHarnessAdapter


class ScannerComponentHarness(ComponentHarnessAdapter, Component):


    class Inventory(Component.Inventory):

        import pyre.inventory
        MacroString = pyre.inventory.str
        InputFile = pyre.inventory.str

        # Basename for all files (may be overridden by specific filename entries).
        fileRoot = MacroString("fileRoot", default="pt1")
        fileRoot.meta['tip'] = "Root pathname for simulation (all filenames derive from this)."
        inputFileRoot = MacroString("inputFileRoot", default="${fileRoot}")
        inputFileRoot.meta['tip'] = "Root input pathname for simulation (all input filenames derive from this)."
        outputFileRoot = MacroString("outputFileRoot", default="${fileRoot}.${rank}")
        outputFileRoot.meta['tip'] = "Root output pathname for simulation (all output filenames derive from this)."
        
        # Optional input files.
        keywordEqualsValueFile = InputFile("keywordEqualsValueFile",default="${inputFileRoot}.keyval")
        keywordEqualsValueFile.meta['tip'] = "Pathname for keyword = value file (overrides default from inputFileRoot)."


    def __init__(self):
        Component.__init__(self, "pl3dscan", "scanner")
        ComponentHarnessAdapter.__init__(self)


    def _configure(self):
        super(ScannerComponentHarness, self)._configure()
        self.fileRoot = self.inventory.fileRoot
        self.inputFileRoot = self.inventory.fileRoot
        self.outputFileRoot = self.inventory.outputFileRoot
        self.keywordEqualsValueFile = self.inventory.keywordEqualsValueFile


    def _init(self):
        super(ScannerComponentHarness, self)._init()
        self.scanner = self.harnessComponent()


    def createComponent(self):
        """create the harnessed component"""

        return Scanner()


    def prepareComponentConfiguration(self, component):
        """prepare the settings for the harnessed component"""

        # Expand macros.

        macros = self.getMacros()
        component.defineMacros(macros)

        # Read traits from the .keyval file.
        
        from pyre.util import expandMacros
        from CodecKeyVal import CodecKeyVal
        from os.path import splitext

        filename = expandMacros(self.keywordEqualsValueFile, macros)
        base, ext = splitext(filename)
        assert ext == ".keyval"
        codec = CodecKeyVal()
        shelf = codec.open(base, 'r')

        if shelf:
            registry = shelf['inventory']
            if registry:
                self.componentRegistry.update(registry)

        return super(ScannerComponentHarness, self).prepareComponentConfiguration(component)


    def getMacros(self):
        from pyre.util import expandMacros

        macros = {
            'fileRoot': self.fileRoot,
            }

        # Expand occurrences of '${fileRoot}'.  (Expansion of other
        # macros is deferred until later.)
        macros['inputFileRoot'] = expandMacros(self.inputFileRoot, macros)
        macros['outputFileRoot'] = expandMacros(self.outputFileRoot, macros)

        return macros


Pylith3d_scan = ScannerComponentHarness


# version
# $Id: Pylith3d_scan.py,v 1.19 2005/06/24 20:22:03 willic3 Exp $

# End of file 
