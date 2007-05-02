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


from cig.cs.petsc import PetscApplication
import os
import PyLithLib
import PyLithMeshLib

from pyre.units.pressure import Pa
from pyre.units.length import m
from pyre.units.time import s


QuadratureOrder = dict(
    Full      = 1,
    Reduced   = 2,
    Selective = 3,
    )


AsciiOutput = dict(
    none = 0,
    echo = 1,
    full = 2,
    )


PlotOutput = dict(
    none   = 0,
    ascii  = 1,
    binary = 2,
    )


AnalysisType = dict(
    dataCheck       = 0,
    stiffnessFactor = 1,
    elasticSolution = 2,
    fullSolution    = 3,
    )


def enumConverter(enum):
    def converter(pylithApp, value):
        return enum[value]
    return converter


def ucdOutputConverter(pylithApp, ucdOutput):
    binIOError = None
    try:
        PyLithLib.try_binio(pylithApp.kucd)
    except RuntimeError, binIOError:
        iucd = 1
    else:
        iucd = 2
    if ucdOutput == "none":
        iucd = 0
    elif ucdOutput == "ascii":
        iucd = 1
    elif ucdOutput == "binary":
        if binIOError is None:
            iucd = 2
        else:
            import journal
            warning = journal.warning("pylith3d")
            warning.line("Forcing 'ucdOutput' to 'ascii'.")
            warning.line("Binary UCD output not supported for this Fortran compiler.")
            warning.log(binIOError)
    return iucd


identity = lambda self, value: value


class PyLithApp(PetscApplication):


    name = "pylith3d"


    #
    # properties
    #

    import pyre.inventory as pyre

    MacroString = pyre.str
    OutputFile = pyre.str
    InputFile = pyre.str

    # Title
    title = pyre.str("title", default="PyLith-0.8 Simulation")
    title.meta['tip'] = "Title for this simulation"
    title.__converter = identity

    # Basename for all files (may be overridden by specific filename entries).
    fileRoot = pyre.str("fileRoot", default="pt1")
    fileRoot.meta['tip'] = "Root pathname for simulation (all filenames derive from this)."
    fileRoot.__converter = identity
    inputFileRoot = pyre.str("inputFileRoot", default="${fileRoot}")
    inputFileRoot.meta['tip'] = "Root input pathname for simulation (all input filenames derive from this)."
    outputFileRoot = pyre.str("outputFileRoot", default="${fileRoot}")
    outputFileRoot.meta['tip'] = "Root output pathname for simulation (all output filenames derive from this)."

    # Output filenames (all are optional).
    ofile = OutputFile("asciiOutputFile",default="${outputFileRoot}.ascii")
    ofile.meta['tip'] = "Pathname for ascii output file (overrides default from outputFileRoot)."
    ofile.__converter = identity

    pfile = OutputFile("plotOutputFile",default="${outputFileRoot}.plot")
    pfile.meta['tip'] = "Pathname for plot output file (overrides default from outputFileRoot)."
    pfile.__converter = identity

    ucdroot = MacroString("ucdOutputRoot",default="${outputFileRoot}")
    ucdroot.meta['tip'] = "Base name for UCD output files (overrides default from outputFileRoot)."
    ucdroot.__converter = identity

    # Required input files.
    coordinateInputFile = InputFile("coordinateInputFile",default="${inputFileRoot}.coord")
    coordinateInputFile.meta['tip'] = "Pathname for coordinate input file (overrides default from inputFileRoot)."
    coordinateInputFile.__converter = identity

    bcfile = InputFile("bcInputFile",default="${inputFileRoot}.bc")
    bcfile.meta['tip'] = "Pathname for boundary condition input file (overrides default from inputFileRoot)."
    bcfile.__converter = identity

    timeStepInputFile = InputFile("timeStepInputFile",default="${inputFileRoot}.time")
    timeStepInputFile.meta['tip'] = "Pathname for time step definitions input file (overrides default from inputFileRoot)."
    timeStepInputFile.__converter = identity

    stfile = InputFile("stateVariableInputFile",default="${inputFileRoot}.statevar")
    stfile.meta['tip'] = "Pathname for file defining which state variables to output (overrides default from inputFileRoot)."
    stfile.__converter = identity

    materialPropertiesInputFile = InputFile("materialPropertiesInputFile",default="${inputFileRoot}.prop")
    materialPropertiesInputFile.meta['tip'] = "Pathname for file defining material properties (overrides default from inputFileRoot)."

    connectivityInputFile = InputFile("connectivityInputFile",default="${inputFileRoot}.connect")
    connectivityInputFile.meta['tip'] = "Pathname for connectivity input file (overrides default from inputFileRoot)."
    connectivityInputFile.__converter = identity

    # This file is only required for time-dependent problems.
    fofile = InputFile("fullOutputInputFile",default="${inputFileRoot}.fuldat")
    fofile.meta['tip'] = "Pathname for file defining when to provide output (overrides default from inputFileRoot)."
    fofile.__converter = identity

    # These files are optional unless generating Green's functions, in which case they are required.
    sampleLocationFile = InputFile("sampleLocationFile",default="${inputFileRoot}.sample")
    sampleLocationFile.meta['tip'] = "Pathname for Green's function sample locations (overrides default from inputFileRoot)."
    sampleLocationFile.__converter = identity

    spfile = InputFile("splitNodeInputFile",default="${inputFileRoot}.split")
    spfile.meta['tip'] = "Pathname for split node input file (overrides default from inputFileRoot)."
    spfile.__converter = identity

    # Optional input files.
    skfile = InputFile("rotationInputFile",default="${inputFileRoot}.skew")
    skfile.meta['tip'] = "Pathname for skew rotations input file (overrides default from inputFileRoot)."
    skfile.__converter = identity

    hfile = InputFile("loadHistoryInputFile",default="${inputFileRoot}.hist")
    hfile.meta['tip'] = "Pathname for file defining load histories (overrides default from inputFileRoot)."
    hfile.__converter = identity

    tractionInputFile = InputFile("tractionInputFile",default="${inputFileRoot}.traction")
    tractionInputFile.meta['tip'] = "Pathname for traction BC input file (overrides default from inputFileRoot)."
    tractionInputFile.__converter = identity

    # Unused input files.
    wfile = InputFile("winklerInputFile",default="${inputFileRoot}.wink")
    wfile.meta['tip'] = "Pathname for Winkler force input file (overrides default from inputFileRoot)."
    wfile.__converter = identity

    materialHistoryInputFile = InputFile("materialHistoryInputFile",default="${inputFileRoot}.mhist")
    materialHistoryInputFile.meta['tip'] = "Pathname for file defining material histories (overrides default from inputFileRoot -- presently unused)."
    materialHistoryInputFile.__converter = identity

    prestressInputFile = InputFile("prestressInputFile",default="${inputFileRoot}.prestr")
    prestressInputFile.meta['tip'] = "Pathname for prestress input file (overrides default from inputFileRoot -- presently unused)."
    prestressInputFile.__converter = identity

    slfile = InputFile("slipperyNodeInputFile",default="${inputFileRoot}.slip")
    slfile.meta['tip'] = "Pathname for slippery node input file (overrides default from inputFileRoot -- presently unused)."
    slfile.__converter = identity

    difile = InputFile("differentialForceInputFile",default="${inputFileRoot}.diff")
    difile.meta['tip'] = "Pathname for file defining slippery node differential forces (overrides default from inputFileRoot -- presently unused)."
    difile.__converter = identity

    wxfile = InputFile("slipperyWinklerInputFile",default="${inputFileRoot}.winkx")
    wxfile.meta['tip'] = "Pathname for file defining slippery node Winkler forces (overrides default from inputFileRoot -- presently unused)."
    wxfile.__converter = identity

    # Output option flags.
    idout = pyre.str("asciiOutput",default="echo")
    idout.validator = pyre.choice(["none","echo","full"])
    idout.meta['tip'] = "Type of ascii output desired (none, echo, full)."
    idout.__converter = enumConverter(AsciiOutput)

    idsk = pyre.str("plotOutput",default="none")
    idsk.validator = pyre.choice(["none","ascii","binary"])
    idsk.meta['tip'] = "Type of plot output desired (none, ascii, binary)."
    idsk.__converter = enumConverter(PlotOutput)

    iucd = pyre.str("ucdOutput",default=None) # default is 'binary' if available; 'ascii' otherwise
    iucd.validator = pyre.choice(["none","ascii","binary"])
    iucd.meta['tip'] = "Type of UCD output desired (none, ascii, binary)."
    iucd.__converter = ucdOutputConverter

    # Additional option flags.
    icode = pyre.str("analysisType",default="fullSolution")
    icode.validator = pyre.choice(["dataCheck","stiffnessFactor",
                                   "elasticSolution","fullSolution"])
    icode.meta['tip'] = "Type of analysis (dataCheck, stiffnessFactor, elasticSolution, fullSolution)."
    icode.__converter = enumConverter(AnalysisType)

    pythonTimestep = pyre.bool("pythonTimestep",default=False)
    pythonTimestep.meta['tip'] = "Whether to use python timestepping loop (enables VTK output for time-dependent solution)."
    pythonTimestep.__converter = identity

    generateGreen = pyre.bool("generateGreen",default=False)
    generateGreen.meta['tip'] = "Whether to generate Green's function results for the specified split node inputs."
    generateGreen.__converter = identity

    idebug = pyre.bool("debuggingOutput",default=False)
    idebug.meta['tip'] = "Whether to produce debugging output."
    idebug.__converter = lambda self, value: int(value)

    ncycle = pyre.int("numberCycles",default=1)
    ncycle.meta['tip'] = "Number of cycles of the given timestep definitions to perform (default=1)."
    ncycle.__converter = identity

    interpolateMesh = pyre.bool("interpolateMesh",default=False)
    interpolateMesh.meta['tip'] = "Create intermediate mesh entities, such as edges and faces."

    partitioner = pyre.str("partitioner",default="chaco")
    partitioner.validator = pyre.choice(["chaco","parmetis"])
    partitioner.meta['tip'] = "Partitioner (chaco, parmetis)."

    # Unused option flags.
    iskopt = pyre.bool("autoRotateSlipperyNodes",default=True)
    iskopt.meta['tip'] = "Whether to performa automatic rotation for slippery nodes (presently unused)."
    iskopt.__converter = lambda self, value: int(value) + 1

    #
    # category 2 parameters formerly placed in *.keyval files
    #

    stol = pyre.dimensional("stressTolerance", default=1.0e-12*Pa)
    stol.__converter = lambda self, value: value / Pa  # convert to Pa
    dtol = pyre.float("minimumStrainPerturbation", default=1.0e-7)
    dtol.__converter = identity
    epert = pyre.float("initialStrainPerturbation", default=1.0e-1)
    epert.__converter = identity

    nprevdflag = pyre.int("usePreviousDisplacementFlag", default=0)
    nprevdflag.__converter = identity

    intord = pyre.str("quadratureOrder", default="Full")
    intord.validator = pyre.choice(["Full", "Reduced", "Selective"])
    intord.__converter = enumConverter(QuadratureOrder)

    ipstrs = pyre.bool("prestressAutoCompute", default=False)
    ipstrs.__converter = lambda self, value: int(value)
    ipauto = pyre.bool("prestressAutoChangeElasticProps", default=False)
    ipauto.__converter = lambda self, value: int(value)
    tpois = pyre.float("prestressAutoComputePoisson", default=0.49)
    tpois.__converter = identity
    tyoungs = pyre.dimensional("prestressAutoComputeYoungs", default=1.0e30*Pa)
    tyoungs.__converter = lambda self, value: value / Pa  # convert to Pa

    gravityX = pyre.dimensional("gravityX", default=0.0*m/(s*s))
    gravityY = pyre.dimensional("gravityY", default=0.0*m/(s*s))
    gravityZ = pyre.dimensional("gravityZ", default=0.0*m/(s*s))

    prestressScaleXx = pyre.float("prestressScaleXx", default=1.0)
    prestressScaleYy = pyre.float("prestressScaleYy", default=1.0)
    prestressScaleZz = pyre.float("prestressScaleZz", default=1.0)
    prestressScaleXy = pyre.float("prestressScaleXy", default=1.0)
    prestressScaleXz = pyre.float("prestressScaleXz", default=1.0)
    prestressScaleYz = pyre.float("prestressScaleYz", default=1.0)

    winklerScaleX = pyre.float("winklerScaleX", default=1.0)
    winklerScaleY = pyre.float("winklerScaleY", default=1.0)
    winklerScaleZ = pyre.float("winklerScaleZ", default=1.0)

    winklerSlipScaleX = pyre.float("winklerSlipScaleX", default=1.0)
    winklerSlipScaleY = pyre.float("winklerSlipScaleY", default=1.0)
    winklerSlipScaleZ = pyre.float("winklerSlipScaleZ", default=1.0)

    kti = pyre.int("f77StandardInput", default=5)
    kti.__converter = identity
    kto = pyre.int("f77StandardOutput", default=6)
    kto.__converter = identity
    kr = pyre.int("f77FileInput", default=10)
    kr.__converter = identity
    kw = pyre.int("f77AsciiOutput", default=11)
    kw.__converter = identity
    kp = pyre.int("f77PlotOutput", default=12)
    kp.__converter = identity
    kucd = pyre.int("f77UcdOutput", default=13)
    kucd.__converter = identity




    # Tell the framework where to find PETSc functions.
    import PyLithLib as petsc


    def _defaults(self):
        super(PyLithApp, self)._defaults()
        
        # Set defaults for PETSc options.  These can be overridden in a .cfg file
        # or from the command line using the "petsc" component.
        self.setPetscDefaults(dict(
            ksp_monitor        = "true",
            ksp_view           = "true",
            ksp_rtol           = "1.0e-9",
            log_summary        = "true",
            partitioner        = "chaco",
            pc_type            = "bjacobi",
            sub_pc_type        = "ilu",
        ))
        
        return
    

    def _validate(self, context):

        super(PyLithApp, self)._validate(context)

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

        self.ofile                       = outputFile(Inventory.ofile, self.idout == "none" and unused or required)
        self.pfile                       = outputFile(Inventory.pfile, self.idsk == "none" and unused or required)
        self.coordinateInputFile         = inputFile(Inventory.coordinateInputFile,         required)
        self.bcfile                      = inputFile(Inventory.bcfile,                      required)
        self.wfile                       = inputFile(Inventory.wfile,                       unused)
        self.skfile                      = inputFile(Inventory.skfile,                      optional)
        self.timeStepInputFile           = inputFile(Inventory.timeStepInputFile,           required)
        self.fofile                      = inputFile(Inventory.fofile, self.icode == "fullSolution" and required or unused)
        self.stfile                      = inputFile(Inventory.stfile,                      required)
        self.hfile                       = inputFile(Inventory.hfile,                       optional)
        self.materialPropertiesInputFile = inputFile(Inventory.materialPropertiesInputFile, required)
        self.materialHistoryInputFile    = inputFile(Inventory.materialHistoryInputFile,    unused)
        self.connectivityInputFile       = inputFile(Inventory.connectivityInputFile,       required)
        self.prestressInputFile          = inputFile(Inventory.prestressInputFile,          unused)
        self.tractionInputFile           = inputFile(Inventory.tractionInputFile,           optional)
        self.spfile                      = inputFile(Inventory.spfile, self.generateGreen and required or optional)
        # Slippery nodes are not yet implemented in PyLith-0.8.
        self.slfile                      = inputFile(Inventory.slfile,                      unused)
        self.difile                      = inputFile(Inventory.difile,                      unused)
        self.wxfile                      = inputFile(Inventory.wxfile,                      unused)
        self.sampleLocationFile          = inputFile(Inventory.sampleLocationFile, self.generateGreen and required or unused)
        # The call to glob() is somewhat crude -- basically, determine
        # if any files might be in the way.
        self.ucdroot                     = macroString(Inventory.ucdroot)

        if False: # broken
            from glob import glob
            ucdFiles = ([self.ucdroot + ".mesh.inp",
                         self.ucdroot + ".gmesh.inp",
                         self.ucdroot + ".mesh.time.prest.inp",
                         self.ucdroot + ".gmesh.time.prest.inp"]
                        + glob(self.ucdroot + ".mesh.time.[0-9][0-9][0-9][0-9][0-9].inp")
                        + glob(self.ucdroot + ".gmesh.time.[0-9][0-9][0-9][0-9][0-9].inp"))
            item = Inventory.ucdroot
            for ucdFile in ucdFiles:
                try:
                    stream = os.fdopen(os.open(ucdFile, os.O_WRONLY|os.O_CREAT|os.O_EXCL), "w")
                except (OSError, IOError), error:
                    context.error(error, items=[item])
                    break
                else:
                    stream.close()
                    os.remove(ucdFile)

        return


    def _init(self):
        super(PyLithApp, self)._init()
        
        pylith = PyLithLib.PyLith()
        self.pylith = pylith

        # Copy my traits to native code.
        for trait in self.properties():
            value = self.getTraitValue(trait.name)
            converter = getattr(trait, '_PyLithApp__converter', None)
            if converter:
                value = converter(self, value)
                setattr(pylith, trait.attr, value)

        # Copy traits which have a many-to-one mapping.
        
        pylith.grav = (
            self.gravityX / (m/(s*s)),
            self.gravityY / (m/(s*s)),
            self.gravityZ / (m/(s*s)),
            )
        
        pylith.prscal = (
            self.prestressScaleXx,
            self.prestressScaleYy,
            self.prestressScaleZz,
            self.prestressScaleXy,
            self.prestressScaleXz,
            self.prestressScaleYz,
            )

        pylith.wscal = (
            self.winklerScaleX,
            self.winklerScaleY,
            self.winklerScaleZ,
            )

        pylith.wxscal = (
            self.winklerSlipScaleX,
            self.winklerSlipScaleY,
            self.winklerSlipScaleZ,
            )

        self.readMaterialProperties()

        return


    def readMaterialProperties(self):
        from Materials import Materials
        
        matinfo = Materials()
        pylith = self.pylith
        
        pylith.numat = matinfo.readprop(self.materialPropertiesInputFile)
        pylith.prop = matinfo.propertyList
        pylith.infmat = matinfo.materialModel
        
        return


    def main(self, *args, **kwds):
        
        timer = False
        if timer:
            from time import clock as now
            start = now()

        from mpi import MPI_Comm_rank, MPI_COMM_WORLD
        self.rank = MPI_Comm_rank(MPI_COMM_WORLD)

        points = None
        if self.generateGreen:
            points      = self.readSamplePoints(self.macroString(self.metainventory.sampleLocationFile))

        self.pylith.mesh = PyLithMeshLib.Mesh(self.macroString(self.metainventory.inputFileRoot),
                                              self.macroString(self.metainventory.bcfile),
                                              self.interpolateMesh,
                                              self.partitioner)

        self.initialize()
        
        self.pylith.run(points)

        if timer:
            finish = now()
            usertime = finish - start
            print "Total user time:  %g" % usertime
        self.pylith.mesh = None
        return


    def readSamplePoints(self, filename):
        '''Read in the sampling locations
        - One point per line, three values per line (x,y,z)
        - Returns a Numeric array'''
        import Numeric
        f = file(filename)
        points = []
        for line in f.readlines():
            points.append([float(v) for v in line.strip().split(' ')])
        f.close()
        return Numeric.array(points)


    def initialize(self):

        inputFile = lambda item, category: self.inputFile(item, category, None)
        outputFile = lambda item, category:  self.outputFile(item, category, None)
        macroString = self.macroString
        Inventory = self.metainventory
        optional = self.IOFileCategory(True,   0,      "optional")
        required = self.IOFileCategory(True,   1,       None)

        # Get base file names
        self.ofile                       = outputFile(Inventory.ofile,                      optional)
        self.pfile                       = outputFile(Inventory.pfile,                      optional)
        self.ucdroot                     = macroString(Inventory.ucdroot)
        self.coordinateInputFile         = inputFile(Inventory.coordinateInputFile,         required)
        self.connectivityInputFile       = inputFile(Inventory.connectivityInputFile,       required)
        self.bcfile                      = inputFile(Inventory.bcfile,                      required)
        self.spfile                      = inputFile(Inventory.spfile,                      optional)
        self.tractionInputFile           = inputFile(Inventory.tractionInputFile,           optional)

        # Create filenames for each process
        for attr in ['ofile',
                     'pfile',
                     'ucdroot',
                     'coordinateInputFile',
                     'connectivityInputFile',
                     'bcfile',
                     'spfile',
                     'tractionInputFile']:
            filename = getattr(self, attr)
            s = filename.split('.')
            sieveFilename = ".".join(s[0:1] + [str(self.rank)] + s[1:])
            setattr(self.pylith, attr, sieveFilename)

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
        builtins = {}
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
        if True:
            # This function is out-of-order.
            value = self.macroString(item)
            return value
        value, stream = self.ioFileStream(item, os.O_WRONLY|os.O_CREAT|os.O_EXCL, "w", category, context)
        if stream is not None:
            stream.close()
            os.remove(value)
        return value


# end of file 
