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


class Pylith3d_scan(Component):


    def __init__(self):
        Component.__init__(self, "pl3dscan", "scanner")

        import journal
        self.trace = journal.debug("pylith3d.trace")

        self.trace.log("Hello from pl3dscan.__init__!")
        
        self.rank = 0
        
        return


    def _validate(self, context):

        super(Pylith3d_scan, self)._validate(context)

        #
        # Open input files.  Log I/O errors.
        #
        
        inputFile = lambda item, category: self.inputFile(item, category, context)
        outputFile = lambda item, category:  self.outputFile(item, category, context)
        macroString = self.macroString

        #                              open?   fatal?  label
        optional = self.IOFileCategory(True,   False,  "optional")
        unused   = self.IOFileCategory(False,  False,  "unused")
        required = self.IOFileCategory(True,   True,    None)
        
        Inventory = self.metainventory

        analysisType = self.inventory.analysisType

        self.asciiOutputFile             = outputFile(Inventory.asciiOutputFile,            optional)
        self.plotOutputFile              = outputFile(Inventory.plotOutputFile,             optional)
        self.coordinateInputFile         = inputFile(Inventory.coordinateInputFile,         required)
        self.bcInputFile                 = inputFile(Inventory.bcInputFile,                 required)
        self.winklerInputFile            = inputFile(Inventory.winklerInputFile,            unused)
        self.rotationInputFile           = inputFile(Inventory.rotationInputFile,           optional)
        self.timeStepInputFile           = inputFile(Inventory.timeStepInputFile,           required)
        self.fullOutputInputFile         = inputFile(Inventory.fullOutputInputFile, analysisType == "fullSolution" and required or unused)
        self.stateVariableInputFile      = inputFile(Inventory.stateVariableInputFile,      required)
        self.loadHistoryInputFile        = inputFile(Inventory.loadHistoryInputFile,        optional)
        self.materialPropertiesInputFile = inputFile(Inventory.materialPropertiesInputFile, required)
        self.materialHistoryInputFile    = inputFile(Inventory.materialHistoryInputFile,    unused)
        self.connectivityInputFile       = inputFile(Inventory.connectivityInputFile,       required)
        self.prestressInputFile          = inputFile(Inventory.prestressInputFile,          unused)
        self.tractionInputFile           = inputFile(Inventory.tractionInputFile,           optional)
        self.splitNodeInputFile          = inputFile(Inventory.splitNodeInputFile,          optional)
        # Slippery nodes are not yet implemented in PyLith-0.8.
        self.slipperyNodeInputFile       = inputFile(Inventory.slipperyNodeInputFile,       unused)
        self.differentialForceInputFile  = inputFile(Inventory.differentialForceInputFile,  unused)
        self.slipperyWinklerInputFile    = inputFile(Inventory.slipperyWinklerInputFile,    unused)
        self.sampleLocationFile          = inputFile(Inventory.sampleLocationFile,          optional)

        # The call to glob() is somewhat crude -- basically, determine
        # if any files might be in the way.
        self.ucdOutputRoot               = macroString(Inventory.ucdOutputRoot)

        if False: # broken
            from glob import glob
            ucdFiles = ([self.ucdOutputRoot + ".mesh.inp",
                         self.ucdOutputRoot + ".gmesh.inp",
                         self.ucdOutputRoot + ".mesh.time.prest.inp",
                         self.ucdOutputRoot + ".gmesh.time.prest.inp"]
                        + glob(self.ucdOutputRoot + ".mesh.time.[0-9][0-9][0-9][0-9][0-9].inp")
                        + glob(self.ucdOutputRoot + ".gmesh.time.[0-9][0-9][0-9][0-9][0-9].inp"))
            item = Inventory.ucdOutputRoot
            for ucdFile in ucdFiles:
                try:
                    stream = os.fdopen(os.open(ucdFile, os.O_WRONLY|os.O_CREAT|os.O_EXCL), "w")
                except (OSError, IOError), error:
                    context.error(error, items=[item])
                    break
                else:
                    stream.close()
                    os.remove(ucdFile)

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
        
        self.rank = MPI_Comm_rank(MPI_COMM_WORLD)
        inputFile = lambda item, category: self.inputFile(item, category, None)
        outputFile = lambda item, category:  self.outputFile(item, category, None)
        macroString = self.macroString
        Inventory = self.metainventory
        optional = self.IOFileCategory(True,   0,      "optional")
        required = self.IOFileCategory(True,   1,       None)

        self.trace.log("Hello from pl3dscan.initialize (begin)!")
        
        self.trace.log("Scanning ascii files to determine dimensions:")

        # Get base file names
        self.asciiOutputFile             = outputFile(Inventory.asciiOutputFile,            optional)
        self.plotOutputFile              = outputFile(Inventory.plotOutputFile,             optional)
        self.ucdOutputRoot               = macroString(Inventory.ucdOutputRoot)
        self.coordinateInputFile         = inputFile(Inventory.coordinateInputFile,         required)
        self.connectivityInputFile       = inputFile(Inventory.connectivityInputFile,       required)
        self.bcInputFile                 = inputFile(Inventory.bcInputFile,                 required)
        self.splitNodeInputFile          = inputFile(Inventory.splitNodeInputFile,          optional)
        self.tractionInputFile           = inputFile(Inventory.tractionInputFile,           optional)

        # Create filenames for each process
        for attr in ['asciiOutputFile',
                     'plotOutputFile',
                     'ucdOutputRoot',
                     'coordinateInputFile',
                     'connectivityInputFile',
                     'bcInputFile',
                     'splitNodeInputFile',
                     'tractionInputFile']:
            filename = getattr(self, attr)
            s = filename.split('.')
            sieveFilename = ".".join(s[0:1] + [str(self.rank)] + s[1:])
            setattr(self, attr + 'Sieve', sieveFilename)

        uparser = pyre.units.parser()
        matinfo = Materials()


        # Define information needed from other functions:
        f77FileInput = self.f77FileInput
        prestressAutoCompute = self.prestressAutoCompute
        prestressAutoChangeElasticProps = self.prestressAutoChangeElasticProps
        quadratureOrder = self.quadratureOrder

        # Initialization of all parameters
	# Memory size variable to keep approximate track of all
	# allocated memory.  This does not include python variables and
	# lists.
	self.memorySize = 0L
	self.intSize = 4L
	self.doubleSize = 8L
        # Parameters that are invariant for this geometry type
        self.geometryType = ""
        self.geometryTypeInt = 0
        self.numberSpaceDimensions = 0
        self.numberDegreesFreedom = 0
        # Note:  eventually the variable below should disappear, and the
        # total size of the state variable array for each material model
        # should be used instead.  This means that all state variable
        # bookkeeping should be done within the material model routines.
        self.stateVariableDimension = 0
        self.materialMatrixDimension = 0
        self.numberSkewDimensions = 0
        self.numberSlipDimensions = 0
        self.numberSlipNeighbors = 0
        self.listIddmat = [0]

        # Invariant parameters related to element type
        self.numberElementTypes = 0
        self.numberElementTypesBase = 0
        self.numberElementNodesBase = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.pointerToListArrayNumberElementNodesBase = None

        # Invariant parameters related to material model
        self.maxMaterialModels = 0
        self.maxStateVariables = 0
        self.maxState0Variables = 0
        self.pointerToMaterialModelInfo = None

        # Parameters derived from values in the inventory or the
        # category 2 parameters above.
        self.analysisTypeInt = 0
        self.pythonTimestep = 0
        self.prestressAutoComputeInt = 0
        self.prestressAutoChangeElasticPropsInt = 0
        self.pointerToSh = None
        self.pointerToShj = None
        self.pointerToGauss = None
        self.pointerToSh2d = None
        self.pointerToGauss2d = None

        # Parameters derived from the number of entries in a file

        self.numberNodes = 0
        self.coordinateUnits = "coordinateUnitsInitial12345678"
        self.coordinateScaleFactor = 0.0

        self.numberBcEntries = 0
        self.displacementUnits = "displacementUnitsInitial123456"
        self.displacementScaleFactor = 0.0
        self.velocityUnits = "velocityUnitsInitial1234567890"
        self.velocityScaleFactor = 0.0
        self.forceUnits = "forceUnitsInitial1234567890123"
        self.forceScaleFactor = 0.0

	self.winklerInfo = [0, 0]
        self.numberWinklerEntries = 0
        self.numberWinklerForces = 0

        self.numberRotationEntries = 0
        self.rotationUnits = "rotationUnitsInitial1234567890"
        self.rotationScaleFactor = 0.0

        self.timeStepInfo = [0, 0]
        self.numberTimeStepGroups = 0
        self.totalNumberTimeSteps = 0
        self.timeUnits = "timeUnitsInitial12345678901234"
        self.timeScaleFactor = 0.0

        self.numberFullOutputs = 0

        self.numberLoadHistories = 0

        self.numberMaterials = 0
        self.propertyListSize = 0
        self.propertyList = [0]
        self.pointerToListArrayPropertyList = None
        self.propertyListIndex = [0]
        self.pointerToListArrayPropertyListIndex = None
        self.materialModel = [0]
        self.pointerToListArrayMaterialModel = None

        self.volumeElementDimens = [0, 0, 0]
        self.numberVolumeElements = 0
        self.volumeElementType = 0
        self.numberVolumeElementFamilies = 0
        self.maxNumberVolumeElementFamilies = 0
        self.numberAllowedVolumeElementTypes = 0
        self.pointerToVolumeElementFamilyList = None

        self.numberPrestressEntries = 0

        self.numberTractionBc = 0
        self.tractionBcUnits = "tractionBcUnitsInitial12345678"
        self.tractionBcScaleFactor = 0.0
        self.tractionFlag = 0

        self.numberSplitNodeEntries = 0

        self.numberSlipperyNodeEntries = 0
        self.numberDifferentialForceEntries = 0
	self.slipperyWinklerInfo = [0, 0]
        self.numberSlipperyWinklerEntries = 0
        self.numberSlipperyWinklerForces = 0

        # This is a test version where the geometry type is automatically
        # specified by using Pylith3d.  The geometry type is only used for
        # f77 routines and not in pyre. An integer value is also defined
        # for use in f77 routines.
        # Define some integer values that are derived from string variables.

        # Parameters that are invariant for this geometry type
        self.geometryType = "3D"
        self.geometryTypeInt = 4
        self.numberSpaceDimensions = 3
        self.numberDegreesFreedom = 3
        self.stateVariableDimension = 6
        self.materialMatrixDimension = 21
        self.numberSkewDimensions = 2
        self.numberSlipDimensions = 5
        self.numberSlipNeighbors = 4
        # self.listIddmat = [
        #     1, 2, 3, 4, 5, 6,
        #     2, 7, 8, 9,10,11,
        #     3, 8,12,13,14,15,
        #     4, 9,13,16,17,18,
        #     5,10,14,17,19,20,
        #     6,11,15,18,20,21]
        # Changed this to correspond to BLAS packed symmetric matrix format.
        self.listIddmat = [
             1, 2, 4, 7,11,16,
             2, 3, 5, 8,12,17,
             4, 5, 6, 9,13,18,
             7, 8, 9,10,14,19,
            11,12,13,14,15,20,
            16,17,18,19,20,21]

        # Invariant parameters related to element type
        self.maxElementNodes = 20
        self.maxGaussPoints = 27
        self.maxElementEquations = self.numberDegreesFreedom*self.maxElementNodes
        self.numberElementTypes = 62
        self.numberElementTypesBase = 10
        self.numberElementNodesBase = [8, 7, 6, 5, 4, 20, 18, 15, 13, 10]
        self.pointerToListArrayNumberElementNodesBase = pylith3d.intListToArray(
            self.numberElementNodesBase)
	self.memorySize += self.numberElementTypesBase*self.intSize
        self.maxElementNodes2d = 4
        self.maxGaussPoints2d = 4
        self.numberElementTypes2d = 2
        self.numberElementTypesBase2d = 2
        self.numberElementNodesBase2d = [4, 3]
        self.pointerToListArrayNumberElementNodesBase2d = pylith3d.intListToArray(
            self.numberElementNodesBase2d)
	self.memorySize += self.numberElementTypesBase2d*self.intSize

        # Invariant parameters related to material model
        self.maxMaterialModels = 20
        self.maxStateVariables = 30
        self.maxState0Variables = 6
        self.pointerToMaterialModelInfo = pylith3d.allocateInt(
            6*self.maxMaterialModels)
	self.memorySize += 6*self.maxMaterialModels*self.intSize

        pylith3d.matmod_def(
            self.pointerToMaterialModelInfo)

        # Parameters derived from values in the inventory or the
        # category 2 parameters above.
        analysisType = self.inventory.analysisType
        analysisTypeMap = {
            "dataCheck":       0,
            "stiffnessFactor": 1,
            "elasticSolution": 2,
            "fullSolution":    3,
            }
        self.analysisTypeInt = analysisTypeMap[analysisType]

        if prestressAutoCompute:
            self.prestressAutoComputeInt = 1
        else:
            self.prestressAutoComputeInt = 0

        if prestressAutoChangeElasticProps:
            self.prestressAutoChangeElasticPropsInt = 1
        else:
            self.prestressAutoChangeElasticPropsInt = 0

        # Parameters derived from the number of entries in a file.
        self.numberNodes = pylith3d.scan_coords(
            f77FileInput,
            self.coordinateUnits,
            self.coordinateInputFileSieve)

        self.coordinateScaleString = \
                                    uparser.parse(string.strip(self.coordinateUnits))
        self.coordinateScaleFactor = self.coordinateScaleString.value

        self.numberBcEntries = pylith3d.scan_bc(
            f77FileInput,
            self.displacementUnits,
            self.velocityUnits,
            self.forceUnits,
            self.bcInputFileSieve)

        if self.numberBcEntries > 0:
            self.displacementScaleString = \
                                          uparser.parse(string.strip(self.displacementUnits))
            self.displacementScaleFactor = self.displacementScaleString.value
            self.velocityScaleString = \
                                      uparser.parse(string.strip(self.velocityUnits))
            self.velocityScaleFactor = self.velocityScaleString.value
            self.forceScaleString = \
                                   uparser.parse(string.strip(self.forceUnits))
            self.forceScaleFactor = self.forceScaleString.value

        self.winklerInfo = pylith3d.scan_wink(
            f77FileInput,
            self.winklerInputFile)
        self.numberWinklerEntries = self.winklerInfo[0]
        self.numberWinklerForces = self.winklerInfo[1]

        self.numberRotationEntries = pylith3d.scan_skew(
            f77FileInput,
            self.rotationUnits,
            self.rotationInputFile)

        if self.numberRotationEntries != 0:
            self.rotationScaleString = \
                                      uparser.parse(string.strip(self.rotationUnits))
            self.rotationScaleFactor = self.rotationScaleString.value

        self.timeStepInfo = pylith3d.scan_timdat(
            f77FileInput,
            self.timeUnits,
            self.timeStepInputFile)
        self.numberTimeStepGroups = self.timeStepInfo[0]
        self.totalNumberTimeSteps = self.timeStepInfo[1]

        self.timeScaleString = \
                              uparser.parse(string.strip(self.timeUnits))
        self.timeScaleFactor = self.timeScaleString.value

        self.numberFullOutputs = pylith3d.scan_fuldat(
            self.analysisTypeInt,
            self.totalNumberTimeSteps,
            f77FileInput,
            self.fullOutputInputFile)

        self.numberLoadHistories = pylith3d.scan_hist(
            f77FileInput,
            self.loadHistoryInputFile)

        self.numberMaterials = matinfo.readprop(self.materialPropertiesInputFile)

        self.propertyList = matinfo.propertyList
        self.propertyListIndex = matinfo.propertyIndex
        self.materialModel = matinfo.materialModel
        self.propertyListSize = len(self.propertyList)
        self.pointerToListArrayPropertyList = pylith3d.doubleListToArray(
            self.propertyList)
        self.memorySize += self.propertyListSize*self.doubleSize
        self.pointerToListArrayPropertyListIndex = pylith3d.intListToArray(
            self.propertyListIndex)
        self.memorySize += self.numberMaterials*self.intSize
        self.pointerToListArrayMaterialModel = pylith3d.intListToArray(
            self.materialModel)
        self.memorySize += self.numberMaterials*self.intSize

        # At present, we assume that the number of element families is equal to
        # the number of material types used, since only one volume element type at a
        # time is allowed.
        self.numberAllowedVolumeElementTypes = 1
        self.maxNumberVolumeElementFamilies = self.numberAllowedVolumeElementTypes* \
                                               self.numberMaterials

        self.pointerToVolumeElementFamilyList = pylith3d.allocateInt(
            3*self.maxNumberVolumeElementFamilies)
        self.memorySize += 3*self.maxNumberVolumeElementFamilies*self.intSize

        self.volumeElementDimens = pylith3d.scan_connect(
            self.pointerToListArrayNumberElementNodesBase,
            self.pointerToMaterialModelInfo,
            self.pointerToListArrayMaterialModel,
            self.pointerToVolumeElementFamilyList,
            self.maxNumberVolumeElementFamilies,
	    self.numberMaterials,
            f77FileInput,
            self.connectivityInputFileSieve)

        self.numberVolumeElements = self.volumeElementDimens[0]
        self.numberVolumeElementFamilies = self.volumeElementDimens[1]
        self.volumeElementType = self.volumeElementDimens[2]

        self.pointerToListArrayMaterialModel = None
        self.pointerToListArrayPropertyListIndex = None
        self.memorySize -= 2*self.numberMaterials*self.intSize

        # self.numberPrestressEntries = pylith3d.scan_prestr(
        #     self.stateVariableDimension,
        #     self.numberPrestressGaussPoints,
        #     self.numberElements,
        #     self.prestressAutoComputeInt,
        #     f77FileInput,
        #     self.prestressInputFile)

        self.numberTractionBc = pylith3d.scan_tractions(
            self.maxElementNodes2d,
            f77FileInput,
            self.tractionBcUnits,
            self.tractionInputFileSieve)

        if self.numberTractionBc != 0:
            self.tractionBcScaleString = \
                                        uparser.parse(string.strip(self.tractionBcUnits))
            self.tractionBcScaleFactor = self.tractionBcScaleString.value
            self.tractionFlag = 1

        self.numberSplitNodeEntries = pylith3d.scan_split(
            f77FileInput,
            self.splitNodeInputFileSieve)

        self.numberSlipperyNodeEntries = pylith3d.scan_slip(
            f77FileInput,
            self.slipperyNodeInputFile)

        self.numberDifferentialForceEntries = pylith3d.scan_diff(
            self.numberSlipperyNodeEntries,
            f77FileInput,
            self.differentialForceInputFile)

        self.slipperyWinklerInfo = pylith3d.scan_winkx(
            self.numberSlipperyNodeEntries,
            f77FileInput,
            self.slipperyWinklerInputFile)
        self.numberSlipperyWinklerEntries = self.slipperyWinklerInfo[0]
        self.numberSlipperyWinklerForces = self.slipperyWinklerInfo[1]

        self.trace.log("Hello from pl3dscan.initialize (end)!")

        return


    class IOFileCategory(object):
        def __init__(self, tryOpen, isFatal, label):
            self.tryOpen = tryOpen
            self.isFatal = isFatal
            self.label = label
    
    def macroString(self, item):
        from pyre.util import expandMacros
        class InventoryAdapter(object):
            def __init__(self, inventory, builtins):
                self.inventory = inventory
                self.builtins = builtins
            def __getitem__(self, key):
                builtin = self.builtins.get(key)
                if builtin is None:
                    return expandMacros(str(self.inventory.getTraitValue(key)), self)
                return builtin
        builtins = {
            'rank': str(self.rank),
            }
        return expandMacros(item.value, InventoryAdapter(self.inventory, builtins))

    def ioFileStream(self, item, flags, mode, category, context):
        value = self.macroString(item)
        stream = None
        if category.tryOpen:
            try:
                stream = os.fdopen(os.open(value, flags), mode)
            except (OSError, IOError), error:
                if context is None:
                    if category.isFatal:
                        raise
                elif category.isFatal:
                    context.error(error, items=[item])
                else:
                    pass # warning?
        return value, stream

    def inputFile(self, item, category, context):
        value, stream = self.ioFileStream(item, os.O_RDONLY, "r", category, context)
        if stream is not None:
            stream.close()
        return value
    
    def inputFileStream(self, item, category, context):
        return self.ioFileStream(item, os.O_RDONLY, "r", category, context)[1]
    
    def outputFile(self, item, category, context):
        value, stream = self.ioFileStream(item, os.O_WRONLY|os.O_CREAT|os.O_EXCL, "w", category, context)
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

        # Basename for all files (may be overridden by specific filename entries).
        fileRoot = pyre.inventory.str("fileRoot", default="pt1")
        fileRoot.meta['tip'] = "Root pathname for simulation (all filenames derive from this)."
        inputFileRoot = pyre.inventory.str("inputFileRoot", default="${fileRoot}")
        inputFileRoot.meta['tip'] = "Root input pathname for simulation (all input filenames derive from this)."
        outputFileRoot = pyre.inventory.str("outputFileRoot", default="${fileRoot}")
        outputFileRoot.meta['tip'] = "Root output pathname for simulation (all output filenames derive from this)."
        
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
        rotationInputFile = InputFile("rotationInputFile",default="${inputFileRoot}.skew")
        rotationInputFile.meta['tip'] = "Pathname for skew rotations input file (overrides default from inputFileRoot)."

        loadHistoryInputFile = InputFile("loadHistoryInputFile",default="${inputFileRoot}.hist")
        loadHistoryInputFile.meta['tip'] = "Pathname for file defining load histories (overrides default from inputFileRoot)."

        sampleLocationFile = InputFile("sampleLocationFile",default="${inputFileRoot}.sample")
        sampleLocationFile.meta['tip'] = "Pathname for Green's function sample locations (overrides default from inputFileRoot)."

        splitNodeInputFile = InputFile("splitNodeInputFile",default="${inputFileRoot}.split")
        splitNodeInputFile.meta['tip'] = "Pathname for split node input file (overrides default from inputFileRoot)."

        tractionInputFile = InputFile("tractionInputFile",default="${inputFileRoot}.traction")
        tractionInputFile.meta['tip'] = "Pathname for traction BC input file (overrides default from inputFileRoot)."

        # Unused input files.
        winklerInputFile = InputFile("winklerInputFile",default="${inputFileRoot}.wink")
        winklerInputFile.meta['tip'] = "Pathname for Winkler force input file (overrides default from inputFileRoot)."

        materialHistoryInputFile = InputFile("materialHistoryInputFile",default="${inputFileRoot}.mhist")
        materialHistoryInputFile.meta['tip'] = "Pathname for file defining material histories (overrides default from inputFileRoot -- presently unused)."

        prestressInputFile = InputFile("prestressInputFile",default="${inputFileRoot}.prestr")
        prestressInputFile.meta['tip'] = "Pathname for prestress input file (overrides default from inputFileRoot -- presently unused)."

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
        # category 2 parameters formerly placed in *.keyval files
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

    def setup(self):

        self.trace.log("Hello from pl3dsetup.initialize (begin)!")
        
        # Initialize and define some integer parameters based on string
        # or logical parameters in python

        self.quadratureOrderInt = 0
        if self.quadratureOrder == "Full":
            self.quadratureOrderInt = 1
        elif self.quadratureOrder == "Reduced":
            self.quadratureOrderInt = 2
        elif self.quadratureOrder == "Selective":
            self.quadratureOrderInt = 3
        else:
            self.quadratureOrderInt = 1

        self.asciiOutputInt = 0
        if self.inventory.asciiOutput == "none":
            self.asciiOutputInt = 0
        elif self.inventory.asciiOutput == "echo":
            self.asciiOutputInt = 1
        else:
            self.asciiOutputInt = 2
            
        self.plotOutputInt = 0
        if self.inventory.plotOutput == "none":
            self.plotOutputInt = 0
        elif self.inventory.plotOutput == "ascii":
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
        if self.inventory.ucdOutput == "none":
            self.ucdOutputInt = 0
        elif self.inventory.ucdOutput == "ascii":
            self.ucdOutputInt = 1
        elif self.inventory.ucdOutput == "binary":
            if binIOError is None:
                self.ucdOutputInt = 2
            else:
                import journal
                warning = journal.warning("pylith3d")
                warning.line("Forcing 'ucdOutput' to 'ascii'.")
                warning.line("Binary UCD output not supported for this Fortran compiler.")
                warning.log(binIOError)
            
        self.debuggingOutputInt = 0
        if self.inventory.debuggingOutput:
            self.debuggingOutputInt = 1
        else:
            self.debuggingOutputInt = 0

        self.autoRotateSlipperyNodesInt = 0
        if self.inventory.autoRotateSlipperyNodes:
            self.autoRotateSlipperyNodesInt = 2
        else:
            self.autoRotateSlipperyNodesInt = 1

        # Get some parameters from the inventory list.
        self.title = self.inventory.title
        self.numberCycles = self.inventory.numberCycles

        self.trace.log("Hello from pl3dsetup.initialize (end)!")

        return


    def read(self):

        # This function reads all input and performs some memory allocation.

        from ElementTypeDef import ElementTypeDef
        import pylith3d

        self.trace.log("Hello from pl3dsetup.read (begin)!")
        
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

        # Surface element type information
        self.elementTypeInfo2d = [0, 0, 0, 0]
        self.pointerToListArrayElementTypeInfo2d = None
        self.pointerToSh2d = None
        self.pointerToGauss2d = None
        self.numberSurfaceElementNodes = 0

        # Traction BC
        self.pointerToTractionverts = None
        self.pointerToTractionvals = None

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

        # Set up global integration info.
        eltype.getdef(
            self.volumeElementType,
            self.quadratureOrderInt,
            self.numberSpaceDimensions,
            self.numberDegreesFreedom)

        self.elementTypeInfo = eltype.elementTypeInfo
        self.elementTypeInfo2d = eltype.elementTypeInfo2d
        self.pointerToSh = eltype.pointerToSh
        self.pointerToSh2d = eltype.pointerToSh2d
        self.pointerToShj = eltype.pointerToShj
        self.pointerToGauss = eltype.pointerToGauss
        self.pointerToGauss2d = eltype.pointerToGauss2d
        self.numberVolumeElementNodes = eltype.numberVolumeElementNodes
        self.numberVolumeElementGaussPoints = eltype.numberVolumeElementGaussPoints
        self.numberVolumeElementEquations = eltype.numberVolumeElementEquations
        self.numberSurfaceElementNodes = eltype.numberSurfaceElementNodes
	self.connectivitySize = self.numberVolumeElements*self.numberVolumeElementNodes
        self.pointerToListArrayElementTypeInfo = pylith3d.intListToArray(
            self.elementTypeInfo)
        self.pointerToListArrayElementTypeInfo2d = pylith3d.intListToArray(
            self.elementTypeInfo2d)
	self.memorySize += 8*self.intSize

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

        # Read traction BC
        self.pointerToTractionverts = pylith3d.allocateInt(
            self.numberSurfaceElementNodes*self.numberTractionBc)
        self.pointerToTractionvals = pylith3d.allocateDouble(
            self.numberDegreesFreedom*self.numberTractionBc)

        pylith3d.read_tractions(
            self.pointerToTractionverts,
            self.pointerToTractionvals,
            self.tractionBcScaleFactor,
            self.numberTractionBc,
            self.numberSurfaceElementNodes,
            self.f77FileInput,
            self.tractionInputFile)

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

        self.trace.log("Hello from pl3dsetup.read (end)!")

        return

    def numberequations(self):

        # This functions numbers equations based on BC and slippery node info.

        import pylith3d

        self.trace.log("Hello from pl3dsetup.numberequations (begin)!")
        
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

        self.trace.log("Hello from pl3dsetup.numberequations (end)!")
            
        return
            
        
    def sortmesh(self):

        # This function sorts elements into families and sorts all other items that are
        # affected by this.

        import pylith3d

        self.trace.log("Hello from pl3dsetup.sortmesh (begin)!")
        
        print "Renumbering elements, split nodes, and slippery nodes:"

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
        # in array Iens, and the index array for the new ordering is
        # Indxiel.  The index array for the original ordering is Ielindx.
        # The original element node array (Ien) and the associated
        # material type array (Mat) may be deallocated after sorting.
        self.pointerToIens = pylith3d.allocateInt(
            self.numberVolumeElementNodes*self.numberVolumeElements)
	self.memorySize += self.numberVolumeElementNodes* \
                           self.numberVolumeElements*self.intSize
        self.pointerToIvfamily = pylith3d.allocateInt(
            6*self.numberVolumeElementFamilies)
        self.memorySize += 5*self.numberVolumeElementFamilies*self.intSize

        self.pointerToIvftmp = pylith3d.allocateInt(
            self.numberVolumeElementFamilies)
        self.memorySize += self.numberVolumeElementFamilies*self.intSize

        self.pointerToIndxiel = pylith3d.allocateInt(
            self.numberVolumeElements)
	self.memorySize += self.numberVolumeElements*self.intSize

        self.pointerToIelindx = pylith3d.allocateInt(
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
            self.pointerToIelindx,
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
            
        self.trace.log("Hello from pl3dsetup.sortmesh (end)!")

        return

        
    def sparsesetup(self, mesh):

        # This function sets up sparse matrix and associated storage.

        import pylith3d

        self.trace.log("Hello from pl3dsetup.sparsesetup (begin)!")
        
        print "Setting up sparse matrix storage:"
        
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
        
        self.trace.log("Hello from pl3dsetup.sparsesetup (end)!")

        return
        
    def allocateremaining(self):

        # This function allocates all remaining arrays that are needed for computations.
        
        import pylith3d

        self.trace.log("Hello from pl3dsetup.allocateremaining (begin)!")
        
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

        self.trace.log("Hello from pl3dsetup.allocateremaining (end)!")

        return


    def meshwrite(self):

        # This function outputs mesh information.
        # In the near future, this needs to be broken into classes for
        # Ascii output, plot output, UCD output, etc.

        import pylith3d

        self.trace.log("Hello from pl3dsetup.meshwriteascii (begin)!")
        
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

        # Write traction info
        pylith3d.write_tractions(
            self.pointerToTractionverts,
            self.pointerToTractionvals,
            self.numberTractionBc,
            self.numberSurfaceElementNodes,
            self.f77AsciiOutput,
            self.asciiOutputInt,
            self.asciiOutputFile)
    
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

        self.trace.log("Hello from pl3dsetup.meshwrite (end)!")

        return



# version
# $Id: Pylith3d_scan.py,v 1.19 2005/06/24 20:22:03 willic3 Exp $

# End of file 
