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



# version
# $Id: Pylith3d_scan.py,v 1.19 2005/06/24 20:22:03 willic3 Exp $

# End of file 
