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
import constants
import os
import pylith3d


prestress = False # code for reading prestress input files is presently disabled


class PyLith(PetscApplication):


    name = "pylith3d"


    #
    # properties
    #

    import pyre.inventory as pyre
    from cig.cs.petsc import PetscProperty

    MacroString = pyre.str
    OutputFile = pyre.str
    InputFile = pyre.str

    # declare PETSc options that are of interest to PyLith
    ksp_monitor        = PetscProperty(default="true")
    ksp_view           = PetscProperty(default="true")
    ksp_rtol           = PetscProperty(default="1.0e-9")
    log_summary        = PetscProperty(default="true")
    pc_type            = PetscProperty(default="bjacobi")
    sub_pc_type        = PetscProperty(default="ilu")
    start_in_debugger  = PetscProperty()
    debugger_pause     = PetscProperty()

    # Title
    title = pyre.str("title", default="PyLith-0.8 Simulation")
    title.meta['tip'] = "Title for this simulation"

    # Basename for all files (may be overridden by specific filename entries).
    fileRoot = pyre.str("fileRoot", default="pt1")
    fileRoot.meta['tip'] = "Root pathname for simulation (all filenames derive from this)."
    inputFileRoot = pyre.str("inputFileRoot", default="${fileRoot}")
    inputFileRoot.meta['tip'] = "Root input pathname for simulation (all input filenames derive from this)."
    outputFileRoot = pyre.str("outputFileRoot", default="${fileRoot}")
    outputFileRoot.meta['tip'] = "Root output pathname for simulation (all output filenames derive from this)."

    # Output filenames (all are optional).
    ofile = OutputFile("asciiOutputFile",default="${outputFileRoot}.ascii")
    ofile.meta['tip'] = "Pathname for ascii output file (overrides default from outputFileRoot)."

    pfile = OutputFile("plotOutputFile",default="${outputFileRoot}.plot")
    pfile.meta['tip'] = "Pathname for plot output file (overrides default from outputFileRoot)."

    ucdroot = MacroString("ucdOutputRoot",default="${outputFileRoot}")
    ucdroot.meta['tip'] = "Base name for UCD output files (overrides default from outputFileRoot)."

    # Required input files.
    coordinateInputFile = InputFile("coordinateInputFile",default="${inputFileRoot}.coord")
    coordinateInputFile.meta['tip'] = "Pathname for coordinate input file (overrides default from inputFileRoot)."

    bcfile = InputFile("bcInputFile",default="${inputFileRoot}.bc")
    bcfile.meta['tip'] = "Pathname for boundary condition input file (overrides default from inputFileRoot)."

    timeStepInputFile = InputFile("timeStepInputFile",default="${inputFileRoot}.time")
    timeStepInputFile.meta['tip'] = "Pathname for time step definitions input file (overrides default from inputFileRoot)."

    stfile = InputFile("stateVariableInputFile",default="${inputFileRoot}.statevar")
    stfile.meta['tip'] = "Pathname for file defining which state variables to output (overrides default from inputFileRoot)."

    materialPropertiesInputFile = InputFile("materialPropertiesInputFile",default="${inputFileRoot}.prop")
    materialPropertiesInputFile.meta['tip'] = "Pathname for file defining material properties (overrides default from inputFileRoot)."

    connectivityInputFile = InputFile("connectivityInputFile",default="${inputFileRoot}.connect")
    connectivityInputFile.meta['tip'] = "Pathname for connectivity input file (overrides default from inputFileRoot)."

    # This file is only required for time-dependent problems.
    fofile = InputFile("fullOutputInputFile",default="${inputFileRoot}.fuldat")
    fofile.meta['tip'] = "Pathname for file defining when to provide output (overrides default from inputFileRoot)."

    # These files are optional unless generating Green's functions, in which case they are required.
    sampleLocationFile = InputFile("sampleLocationFile",default="${inputFileRoot}.sample")
    sampleLocationFile.meta['tip'] = "Pathname for Green's function sample locations (overrides default from inputFileRoot)."

    spfile = InputFile("splitNodeInputFile",default="${inputFileRoot}.split")
    spfile.meta['tip'] = "Pathname for split node input file (overrides default from inputFileRoot)."

    # Optional input files.
    skfile = InputFile("rotationInputFile",default="${inputFileRoot}.skew")
    skfile.meta['tip'] = "Pathname for skew rotations input file (overrides default from inputFileRoot)."

    hfile = InputFile("loadHistoryInputFile",default="${inputFileRoot}.hist")
    hfile.meta['tip'] = "Pathname for file defining load histories (overrides default from inputFileRoot)."

    tractionInputFile = InputFile("tractionInputFile",default="${inputFileRoot}.traction")
    tractionInputFile.meta['tip'] = "Pathname for traction BC input file (overrides default from inputFileRoot)."

    # Unused input files.
    wfile = InputFile("winklerInputFile",default="${inputFileRoot}.wink")
    wfile.meta['tip'] = "Pathname for Winkler force input file (overrides default from inputFileRoot)."

    materialHistoryInputFile = InputFile("materialHistoryInputFile",default="${inputFileRoot}.mhist")
    materialHistoryInputFile.meta['tip'] = "Pathname for file defining material histories (overrides default from inputFileRoot -- presently unused)."

    prestressInputFile = InputFile("prestressInputFile",default="${inputFileRoot}.prestr")
    prestressInputFile.meta['tip'] = "Pathname for prestress input file (overrides default from inputFileRoot -- presently unused)."

    slfile = InputFile("slipperyNodeInputFile",default="${inputFileRoot}.slip")
    slfile.meta['tip'] = "Pathname for slippery node input file (overrides default from inputFileRoot -- presently unused)."

    difile = InputFile("differentialForceInputFile",default="${inputFileRoot}.diff")
    difile.meta['tip'] = "Pathname for file defining slippery node differential forces (overrides default from inputFileRoot -- presently unused)."

    wxfile = InputFile("slipperyWinklerInputFile",default="${inputFileRoot}.winkx")
    wxfile.meta['tip'] = "Pathname for file defining slippery node Winkler forces (overrides default from inputFileRoot -- presently unused)."

    # Output option flags.
    asciiOutput = pyre.str("asciiOutput",default="echo")
    asciiOutput.validator = pyre.choice(["none","echo","full"])
    asciiOutput.meta['tip'] = "Type of ascii output desired (none, echo, full)."

    plotOutput = pyre.str("plotOutput",default="none")
    plotOutput.validator = pyre.choice(["none","ascii","binary"])
    plotOutput.meta['tip'] = "Type of plot output desired (none, ascii, binary)."

    ucdOutput = pyre.str("ucdOutput",default=None)
    ucdOutput.validator = pyre.choice(["none","ascii","binary"])
    ucdOutput.meta['tip'] = "Type of UCD output desired (none, ascii, binary)."

    # Additional option flags.
    analysisType = pyre.str("analysisType",default="fullSolution")
    analysisType.validator = pyre.choice(["dataCheck","stiffnessFactor",
                                          "elasticSolution","fullSolution"])
    analysisType.meta['tip'] = "Type of analysis (dataCheck, stiffnessFactor, elasticSolution, fullSolution)."

    pythonTimestep = pyre.bool("pythonTimestep",default=False)
    pythonTimestep.meta['tip'] = "Whether to use python timestepping loop (enables VTK output for time-dependent solution)."

    generateGreen = pyre.bool("generateGreen",default=False)
    generateGreen.meta['tip'] = "Whether to generate Green's function results for the specified split node inputs."

    debuggingOutput = pyre.bool("debuggingOutput",default=False)
    debuggingOutput.meta['tip'] = "Whether to produce debugging output."

    ncycle = pyre.int("numberCycles",default=1)
    ncycle.meta['tip'] = "Number of cycles of the given timestep definitions to perform (default=1)."

    interpolateMesh = pyre.bool("interpolateMesh",default=False)
    interpolateMesh.meta['tip'] = "Create intermediate mesh entities, such as edges and faces."

    partitioner = pyre.str("partitioner",default="chaco")
    partitioner.validator = pyre.choice(["chaco","parmetis"])
    partitioner.meta['tip'] = "Partitioner (chaco, parmetis)."

    # Unused option flags.
    autoRotateSlipperyNodes = pyre.bool("autoRotateSlipperyNodes",default=True)
    autoRotateSlipperyNodes.meta['tip'] = "Whether to performa automatic rotation for slippery nodes (presently unused)."

    #
    # category 2 parameters formerly placed in *.keyval files
    #

    from pyre.units.pressure import Pa
    from pyre.units.length import m
    from pyre.units.time import s

    winklerScaleX = pyre.float("winklerScaleX", default=1.0)
    winklerScaleY = pyre.float("winklerScaleY", default=1.0)
    winklerScaleZ = pyre.float("winklerScaleZ", default=1.0)

    stol = pyre.dimensional("stressTolerance", default=1.0e-12*Pa)
    dtol = pyre.float("minimumStrainPerturbation", default=1.0e-7)
    epert = pyre.float("initialStrainPerturbation", default=1.0e-1)

    nprevdflag = pyre.int("usePreviousDisplacementFlag", default=0)

    quadratureOrder = pyre.str("quadratureOrder", default="Full")
    quadratureOrder.validator = pyre.choice(["Full", "Reduced", "Selective"])

    gravityX = pyre.dimensional("gravityX", default=0.0*m/(s*s))
    gravityY = pyre.dimensional("gravityY", default=0.0*m/(s*s))
    gravityZ = pyre.dimensional("gravityZ", default=0.0*m/(s*s))

    prestressAutoCompute = pyre.bool("prestressAutoCompute", default=False)
    prestressAutoChangeElasticProps = pyre.bool("prestressAutoChangeElasticProps", default=False)
    tpois = pyre.float("prestressAutoComputePoisson", default=0.49)
    tyoungs = pyre.dimensional("prestressAutoComputeYoungs", default=1.0e30*Pa)

    prestressScaleXx = pyre.float("prestressScaleXx", default=1.0)
    prestressScaleYy = pyre.float("prestressScaleYy", default=1.0)
    prestressScaleZz = pyre.float("prestressScaleZz", default=1.0)
    prestressScaleXy = pyre.float("prestressScaleXy", default=1.0)
    prestressScaleXz = pyre.float("prestressScaleXz", default=1.0)
    prestressScaleYz = pyre.float("prestressScaleYz", default=1.0)

    winklerSlipScaleX = pyre.float("winklerSlipScaleX", default=1.0)
    winklerSlipScaleY = pyre.float("winklerSlipScaleY", default=1.0)
    winklerSlipScaleZ = pyre.float("winklerSlipScaleZ", default=1.0)

    f77StandardInput = pyre.int("f77StandardInput", default=5)
    f77StandardOutput = pyre.int("f77StandardOutput", default=6)
    kr = pyre.int("f77FileInput", default=10)
    kw = pyre.int("f77AsciiOutput", default=11)
    kp = pyre.int("f77PlotOutput", default=12)
    kucd = pyre.int("f77UcdOutput", default=13)




    # Tell the framework where to find PETSc functions.
    import pylith3d as petsc


    # Use PETSc-style command line parsing.
    from cig.cs.petsc import PetscCommandlineParser as CommandlineParser


    # hack to recognize old 'pl3dscan.xxx' and 'scanner.xxx' options
    def applyConfiguration(self, context=None):
        # this mimics the standard Pyre order:  <component-name>.xxx overrides <facility-name>.xxx
        for alias in ["scanner", "pl3dscan"]:
            node = self.inventory._priv_registry.extractNode(alias)
            if node:
                node.name = self.name
                self.updateConfiguration(node)
        return super(PyLith, self).applyConfiguration(context)


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


    def outputSampleValues(self, filename, impulseNode, values):
        '''impulse# sample# sample values'''
        # Computing normal to the fault:
        #   Split nodes define the fault
        #   Get all fault faces for a node
        #   Area weighted average of normals
        f = file(filename, 'w')
        for v, values in enumerate(values):
            write(f, '%d %d %g %g %g' % (impulseNode, v, values[0], values[1], values[2]))
        f.close()
        return


    def greenFunction(self, points):
        """
        # Beginning of loop that loops over split node sets, creating
        # an 'impulse' for each one and outputting response values.
        # Below at present is a quasi-C version of the needed code.

        splitField = None
        m = None

        # Need bindings for these
        pylith3d.meshGetSectionPair(mesh, "split", splitField)
        pylith3d.meshGetMesh(self.mesh, m)

        # This is incorrect, but I need something like:
        topology = getTopology(m)
        patch = 0
        eNumbering = pylith3d.getLocalNumbering(topology, patch, ??)
        vNumbering = pylith3d.getLocalNumbering(topology, patch, 0)


        # Need to loop over global nodes
        for node in ??:
            # Need integer and double lists to hold split node info
            faultind = []
            faultvals = []
            indfault = None
            valfault = None
        
            # Need to find elements in the split node 'patch' that contain node.
            numSet = 0
            for elem in ??:
                # Not sure if this does what I want
                if (pylith3d.sieve.baseContains(elem)):
                    # This is totally wrong, but I need to get local element and node numbers,
                    # along with values
                    numSet += 1
                    faultind += [eNumbering.getIndex(elem)]
                    faultind += [vNumbering.getIndex(node)]

                    # Need to look up how to get field values
                    faultvals += [value.x, value.y, value.z]

                    # Create arrays to send to fortran code
                    indfault = intListToArray(faultind)
                    valfault = doubleListToArray(faultvals)

            # Call fortran routine to set specified split values and clear the rest
            pylith3d.setsplit(self.nfault, self.fault, numfn, indfault, valfault, numSet)
            # Solve
            pylith3d.solveElastic()
            values = self.interpolatePoints(points)
            self.outputSampleValues(self.fileRoot+'.output', values)
        return
                    
        
        
                    SectionReal splitField;
                    
                    # Need bindings for this
                    ierr = MeshGetSectionPair(mesh, "split", &splitField);
                    // Loop over split nodes
                    for() {
                    // Loop over elements
                    for() {
                    # Need bindings for this
                    ierr = SectionPairSetFiberDimension(splitField, e, 1);
                    }
                    # Need bindings for this
                    ierr = SectionPairAllocate(splitField);
                    // Loop over elements
                    for() {
                    PetscPair value;

                    value.i = node;
                    value.x = ;
                    value.y = ;
                    value.z = ;
                    # Need bindings for this
                    ierr = SectionPairUpdate(splitField, e, &value);
                    # Major problem right now:  This just updates PETSc/Sieve's copy of splitField.
                    # It does not change the values within PyLith, which have been read from
                    # per-process input files.
                    }
                    // Solve
                    pl3drun.solveElastic()
                    # Need bindings for this
                    ierr = SectionPairClear(splitField);
                    }
                    values = self.interpolatePoints(points)
                    self.outputSampleValues(self.fileRoot+'.output', values)
                    """
        return


    def main(self, *args, **kwds):
    
#        from time import clock as now
#        start = now()

        import journal
        self.trace = journal.debug("pylith3d.trace")
        
        from mpi import MPI_Comm_rank, MPI_COMM_WORLD
        self.rank = MPI_Comm_rank(MPI_COMM_WORLD)

        if self.generateGreen:
            points      = self.readSamplePoints(self.macroString(self.metainventory.sampleLocationFile))

        self.mesh = pylith3d.processMesh(self.macroString(self.metainventory.bcfile),
                                         self.macroString(self.metainventory.inputFileRoot),
                                         self.interpolateMesh,
                                         self.partitioner)

        self.initialize()
        
        self.setup()
        self.read()
        self.numberequations()
        self.sortmesh()
        self.sparsesetup()
        self.allocateremaining()
        self.meshwrite()

        if self.generateGreen:
            self.greenFunction(points)
        else:
            self.runSimulation()
#        finish = now()
#        usertime = finish - start
#        print "Total user time:  %g" % usertime
        return


    def _validate(self, context):

        super(PyLith, self)._validate(context)

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

        self.ofile                       = outputFile(Inventory.ofile,                      optional)
        self.pfile                       = outputFile(Inventory.pfile,                      optional)
        self.coordinateInputFile         = inputFile(Inventory.coordinateInputFile,         required)
        self.bcfile                      = inputFile(Inventory.bcfile,                      required)
        self.wfile                       = inputFile(Inventory.wfile,                       unused)
        self.skfile                      = inputFile(Inventory.skfile,                      optional)
        self.timeStepInputFile           = inputFile(Inventory.timeStepInputFile,           required)
        self.fofile                      = inputFile(Inventory.fofile, self.analysisType == "fullSolution" and required or unused)
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


# derived or automatically-specified quantities (category 3)

    def initialize(self):

        from Materials import Materials
        import pyre.units
        import string
        
        inputFile = lambda item, category: self.inputFile(item, category, None)
        outputFile = lambda item, category:  self.outputFile(item, category, None)
        macroString = self.macroString
        Inventory = self.metainventory
        optional = self.IOFileCategory(True,   0,      "optional")
        required = self.IOFileCategory(True,   1,       None)

        self.trace.log("Hello from pl3dscan.initialize (begin)!")
        
        self.trace.log("Scanning ascii files to determine dimensions:")

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
            setattr(self, attr, sieveFilename)

        uparser = pyre.units.parser()
        matinfo = Materials()


        # poor man's allocation
        coord_units = "coordinateUnitsInitial12345678"
        displacement_units = "displacementUnitsInitial123456"
        velocity_units = "velocityUnitsInitial1234567890"
        force_units = "forceUnitsInitial1234567890123"
        rotation_units = "rotationUnitsInitial1234567890"
        time_units = "timeUnitsInitial12345678901234"
        traction_units = "tractionBcUnitsInitial12345678"

        # This is a test version where the geometry type is automatically
        # specified by using Pylith3d.  The geometry type is only used for
        # f77 routines and not in pyre. An integer value is also defined
        # for use in f77 routines.
        # Define some integer values that are derived from string variables.

        # Invariant parameters related to element type
        self.maxElementEquations = constants.numberDegreesFreedom*constants.maxElementNodes
        self.neni = pylith3d.intListToArray(
            [8, 7, 6, 5, 4, 20, 18, 15, 13, 10])

        # Invariant parameters related to material model
        self.infmatmod = pylith3d.allocateInt(
            6*constants.maxMaterialModels)

        pylith3d.matmod_def(self.infmatmod)

        # Parameters derived from values in the inventory or the
        # category 2 parameters above.
        analysisTypeMap = {
            "dataCheck":       0,
            "stiffnessFactor": 1,
            "elasticSolution": 2,
            "fullSolution":    3,
            }
        self.icode = analysisTypeMap[self.analysisType]

        if self.prestressAutoCompute:
            self.ipstrs = 1
        else:
            self.ipstrs = 0

        if self.prestressAutoChangeElasticProps:
            self.ipauto = 1
        else:
            self.ipauto = 0

        # Parameters derived from the number of entries in a file.
        self.numnp = pylith3d.scan_coords(
            self.kr,
            coord_units,
            self.coordinateInputFile)

        self.cscale = uparser.parse(string.strip(coord_units)).value

        self.numbc = pylith3d.scan_bc(
            self.kr,
            displacement_units,
            velocity_units,
            force_units,
            self.bcfile)

        if self.numbc > 0:
            self.dscale = uparser.parse(string.strip(displacement_units)).value
            self.vscale = uparser.parse(string.strip(velocity_units)).value
            self.fscale = uparser.parse(string.strip(force_units)).value
        else:
            self.dscale = 0.0
            self.vscale = 0.0
            self.fscale = 0.0

        winklerInfo = pylith3d.scan_wink(
            self.kr,
            self.wfile)
        self.nwinke = winklerInfo[0]
        self.nwink = winklerInfo[1]

        self.numrot = pylith3d.scan_skew(
            self.kr,
            rotation_units,
            self.skfile)

        if self.numrot != 0:
            self.runits = uparser.parse(string.strip(rotation_units)).value
        else:
            self.runits = 0.0

        timeStepInfo = pylith3d.scan_timdat(
            self.kr,
            time_units,
            self.timeStepInputFile)
        self.nintg = timeStepInfo[0]
        self.lastep = timeStepInfo[1]

        self.tunits = uparser.parse(string.strip(time_units)).value

        self.icontr = pylith3d.scan_fuldat(
            self.icode,
            self.lastep,
            self.kr,
            self.fofile)

        self.nhist = pylith3d.scan_hist(
            self.kr,
            self.hfile)

        self.numat = matinfo.readprop(self.materialPropertiesInputFile)

        self.materialModel = matinfo.materialModel
        self.prop = pylith3d.doubleListToArray(
            matinfo.propertyList)

        self.scan_connect()

        if prestress:
            self.numberPrestressEntries = pylith3d.scan_prestr(
                constants.stateVariableDimension,
                self.numberPrestressGaussPoints,
                self.numberElements,
                self.ipstrs,
                self.kr,
                self.prestressInputFile)
        else:
            self.numberPrestressEntries = 0

        self.numtractions = pylith3d.scan_tractions(
            constants.nsnodesmax,
            self.kr,
            traction_units,
            self.tractionInputFile)

        if self.numtractions != 0:
            self.tscale = uparser.parse(string.strip(traction_units)).value
        else:
            self.tscale = 0.0

        self.numfn = pylith3d.scan_split(
            self.kr,
            self.spfile)

        self.numslp = pylith3d.scan_slip(
            self.kr,
            self.slfile)

        self.numdif = pylith3d.scan_diff(
            self.numslp,
            self.kr,
            self.difile)

        slipperyWinklerInfo = pylith3d.scan_winkx(
            self.numslp,
            self.kr,
            self.wxfile)
        self.nwinkxe = slipperyWinklerInfo[0]
        self.nwinkx = slipperyWinklerInfo[1]

        self.trace.log("Hello from pl3dscan.initialize (end)!")

        return


    def scan_connect(self):
        infmat = pylith3d.intListToArray(
            self.materialModel)

        # At present, we assume that the number of element families is equal to
        # the number of material types used, since only one volume element type at a
        # time is allowed.
        self.maxvfamilies = constants.numberAllowedVolumeElementTypes* \
                                               self.numat

        self.ivflist = pylith3d.allocateInt(
            3*self.maxvfamilies)

        volumeElementDimens = pylith3d.scan_connect(
            self.neni,
            self.infmatmod,
            infmat,
            self.ivflist,
            self.maxvfamilies,
	    self.numat,
            self.kr,
            self.connectivityInputFile)

        self.numelv = volumeElementDimens[0]
        self.nvfamilies = volumeElementDimens[1]
        self.ietypev = volumeElementDimens[2]

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
        value, stream = self.ioFileStream(item, os.O_WRONLY|os.O_CREAT|os.O_EXCL, "w", category, context)
        if stream is not None:
            stream.close()
            os.remove(value)
        return value


# The main function of this code is to emulate the original functionality of
# input.f in the original version of TECTON.  This code segment controls the
# allocation of memory and the reading of the input file.  Additional functions
# covered by this code include the sparse matrix setup portion, which also does
# some memory allocation.  Additional code sections will call the main elastic
# and time-dependent solution drivers, which are presently f77 subroutines.


    def setup(self):

        self.trace.log("Hello from pl3dsetup.initialize (begin)!")
        
        # Initialize and define some integer parameters based on string
        # or logical parameters in python

        self.intord = 0
        if self.quadratureOrder == "Full":
            self.intord = 1
        elif self.quadratureOrder == "Reduced":
            self.intord = 2
        elif self.quadratureOrder == "Selective":
            self.intord = 3
        else:
            self.intord = 1

        self.idout = 0
        if self.asciiOutput == "none":
            self.idout = 0
        elif self.asciiOutput == "echo":
            self.idout = 1
        else:
            self.idout = 2
            
        self.idsk = 0
        if self.plotOutput == "none":
            self.idsk = 0
        elif self.plotOutput == "ascii":
            self.idsk = 1
        else:
            self.idsk = 2

        binIOError = None
        try:
            pylith3d.try_binio(self.kucd)
        except RuntimeError, binIOError:
            self.iucd = 1
        else:
            self.iucd = 2
        if self.ucdOutput == "none":
            self.iucd = 0
        elif self.ucdOutput == "ascii":
            self.iucd = 1
        elif self.ucdOutput == "binary":
            if binIOError is None:
                self.iucd = 2
            else:
                import journal
                warning = journal.warning("pylith3d")
                warning.line("Forcing 'ucdOutput' to 'ascii'.")
                warning.line("Binary UCD output not supported for this Fortran compiler.")
                warning.log(binIOError)
            
        self.idebug = 0
        if self.debuggingOutput:
            self.idebug = 1
        else:
            self.idebug = 0

        self.iskopt = 0
        if self.autoRotateSlipperyNodes:
            self.iskopt = 2
        else:
            self.iskopt = 1

        self.trace.log("Hello from pl3dsetup.initialize (end)!")

        return


    def read(self):

        # This function reads all input and performs some memory allocation.

        from ElementTypeDef import ElementTypeDef

        self.trace.log("Hello from pl3dsetup.read (begin)!")
        
        print "Reading problem definition and allocating necessary storage:"


        eltype=ElementTypeDef()

        # Make lists that are used as arrays in the f77 function calls below.
        wscal = pylith3d.doubleListToArray(
            [self.winklerScaleX,
             self.winklerScaleY,
             self.winklerScaleZ])

        if prestress:
            self.prscal = pylith3d.doubleListToArray(
                [self.prestressScaleXx,
                 self.prestressScaleYy,
                 self.prestressScaleZz,
                 self.prestressScaleXy,
                 self.prestressScaleXz,
                 self.prestressScaleYz])
        
        wxscal = pylith3d.doubleListToArray(
            [self.winklerSlipScaleX,
             self.winklerSlipScaleY,
             self.winklerSlipScaleZ])

        # Set up global integration info.
        eltype.getdef(
            self.ietypev,
            self.intord)

        self.sh = eltype.sh
        self.sh2d = eltype.sh2d
        self.shj = eltype.shj
        self.gauss = eltype.gauss
        self.gauss2d = eltype.gauss2d
        self.nen = eltype.nen
        self.ngauss = eltype.ngauss
        self.nee = eltype.nee
        self.nsnodes = eltype.nsnodes
	self.connectivitySize = self.numelv*self.nen
        self.infetype = pylith3d.intListToArray(
            eltype.elementTypeInfo)
        self.infetype2d = pylith3d.intListToArray(
            eltype.elementTypeInfo2d)

        # Node-based info (coordinates, displacement arrays, BC, and skew BC).
        self.x = pylith3d.allocateDouble(
            constants.numberSpaceDimensions*self.numnp)
        self.ibond = pylith3d.allocateInt(
            constants.numberDegreesFreedom*self.numnp)
        self.bond = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numnp)
        self.skew = pylith3d.allocateDouble(
            constants.numberSkewDimensions*self.numnp)

        pylith3d.read_coords(
            self.x,
            self.cscale,
            self.numnp,
            self.kr,
            self.coordinateInputFile)

        self.numberConcForces = pylith3d.read_bc(
            self.bond,
            self.dscale,
            self.vscale,
            self.fscale,
            self.ibond,
            self.numnp,
            self.numbc,
            self.kr,
            self.bcfile)

        pylith3d.read_skew(
            self.skew,
            self.runits,
            self.numrot,
            self.numnp,
            self.iskopt,
            self.kr,
            self.skfile)

        # Allocate and read time step, time output, and load history info.
        self.histry = pylith3d.allocateDouble(
            (self.lastep+1)*self.nhist)
        self.maxstp = pylith3d.allocateInt(
            self.nintg)
        self.delt = pylith3d.allocateDouble(
            self.nintg)
        self.alfa = pylith3d.allocateDouble(
            self.nintg)
        self.maxit = pylith3d.allocateInt(
            self.nintg)
        self.ntdinit = pylith3d.allocateInt(
            self.nintg)
        self.lgdef = pylith3d.allocateInt(
            self.nintg)
        self.utol = pylith3d.allocateDouble(
            self.nintg)
        self.ftol = pylith3d.allocateDouble(
            self.nintg)
        self.etol = pylith3d.allocateDouble(
            self.nintg)
        self.itmax = pylith3d.allocateInt(
            self.nintg)
        self.iprint = pylith3d.allocateInt(
            self.icontr)

        # Note that array Times is needed for output, if requested.
        self.times = pylith3d.allocateDouble(
            self.lastep+1)
        self.istatout = pylith3d.allocateInt(
            3*constants.maxStateVariables)
        self.nstatout = pylith3d.allocateInt(3)

        pylith3d.read_timdat(
            self.delt,
            self.alfa,
            self.utol,
            self.ftol,
            self.etol,
            self.times,
            self.tunits,
            self.maxstp,
            self.maxit,
            self.ntdinit,
            self.lgdef,
            self.itmax,
            self.nintg,
            self.lastep,
            self.kr,
            self.timeStepInputFile)

        pylith3d.read_fuldat(
            self.iprint,
            self.icontr,
            self.icode,
            self.ncycle,
            self.lastep,
            self.kr,
            self.fofile)

        pylith3d.read_stateout(
            self.istatout,
            self.nstatout,
            self.kr,
            self.stfile)

        pylith3d.read_hist(
            self.histry,
            self.times,
            self.nhist,
            self.lastep,
            self.kr,
            self.hfile)

        # Allocate and read info on connectivities and prestresses
        self.ien = pylith3d.allocateInt(
            self.nen*self.numelv)
	self.mat = pylith3d.allocateInt(
	    self.numelv)
        if self.numberPrestressEntries != 0 or self.ipstrs != 0:
            self.prestressFlag = 1
        else:
            self.prestressFlag = 0

        pylith3d.read_connect(
            self.ien,
            self.mat,
            self.nen,
            self.numelv,
            self.numnp,
	    self.nvfamilies,
            self.kr,
            self.connectivityInputFile)

        if prestress:
            pylith3d.read_prestr(
                self.stn,
                self.st0,
                self.prscal,
                self.numberStressComponents,
                self.numberGaussPoints,
                self.numberPrestressGaussPoints,
                self.numberElements,
                self.numberPrestressEntries,
                self.ipstrs,
                self.idout,
                self.kr,
                self.kw,
                self.prestressInputFile,
                self.ofile)

        # Read traction BC
        self.tractionverts = pylith3d.allocateInt(
            self.nsnodes*self.numtractions)
        self.tractionvals = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numtractions)

        pylith3d.read_tractions(
            self.tractionverts,
            self.tractionvals,
            self.tscale,
            self.numtractions,
            self.nsnodes,
            self.kr,
            self.tractionInputFile)

        # Read split node info
        self.nfault = pylith3d.allocateInt(
            3*self.numfn)
        self.fault = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numfn)

        self.numflt = pylith3d.read_split(
            self.fault,
            self.nfault,
            self.numfn,
            self.numnp,
            self.numelv,
            self.kr,
            self.spfile)

        # Read slippery node info
        # Note that array Nslip is also required in functions sortmesh and sparsesetup
        # before it can be deallocated.
        self.nslip = pylith3d.allocateInt(
            constants.numberSlipDimensions*self.numslp)
        self.idhist = pylith3d.allocateInt(
            self.numnp)
        self.diforc = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numnp)

        self.numsn = pylith3d.read_slip(
            self.nslip,
            self.numslp,
            self.numnp,
            self.iskopt,
            self.kr,
            self.slfile)

        pylith3d.read_diff(
            self.diforc,
            self.nslip,
            self.idhist,
            self.numslp,
            self.numdif,
            self.numnp,
            self.kr,
            self.difile)
        
        # Read Winkler forces and slippery Winkler forces.
        # All input is finished after this section.
        self.iwinkdef = pylith3d.allocateInt(
            constants.numberDegreesFreedom*self.nwinke)
        self.iwinkid = pylith3d.allocateInt(
            self.nwinke)
        self.winkdef = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.nwinke)

        self.iwinkxdef = pylith3d.allocateInt(
            constants.numberDegreesFreedom*self.nwinkxe)
        self.iwinkxid = pylith3d.allocateInt(
            self.nwinkxe)
        self.winkxdef = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.nwinkxe)

        pylith3d.read_wink(
            self.winkdef,
            wscal,
            self.iwinkdef,
            self.iwinkid,
            self.nwink,
            self.nwinke,
            self.kr,
            self.wfile)

        pylith3d.read_wink(
            self.winkxdef,
            wxscal,
            self.iwinkxdef,
            self.iwinkxid,
            self.nwinkx,
            self.nwinkxe,
            self.kr,
            self.wxfile)

        self.trace.log("Hello from pl3dsetup.read (end)!")

        return

    def numberequations(self):

        # This functions numbers equations based on BC and slippery node info.

        self.trace.log("Hello from pl3dsetup.numberequations (begin)!")
        
        print "Numbering global equations:"

        # Create Idftn array for split nodes.  This can be deallocated after meshwrite function has been called.
        self.idftn = pylith3d.allocateInt(
            self.numflt)

        pylith3d.id_split(
            self.nfault,
            self.idftn,
            self.numnp,
            self.numfn,
            self.numflt)

        # Determine global equations and store equation numbers in Id and Idx.
        self.id = pylith3d.allocateInt(
            constants.numberSpaceDimensions*self.numnp)
        self.idx = pylith3d.allocateInt(
            constants.numberSpaceDimensions*self.numnp)
        self.idslp = pylith3d.allocateInt(
            self.numnp)

        # Number of equations
        self.neq = pylith3d.create_id(
            self.id,
            self.idx,
            self.ibond,
            self.nslip,
            self.idslp,
            self.numslp,
            self.numnp,
            self.numsn)

        self.ipslp = pylith3d.allocateInt(
            constants.numberSlipNeighbors*self.numsn)

        # If there are slippery nodes and the auto-rotation option is selected, find
        # neighboring nodes on the fault so that a best-fit plane can be determined at
        # each node.
        if self.numsn != 0 and self.iskopt ==  2:
            self.nfind()

        # Assign appropriate equation numbers to Iwink array, and compact Wink
        # array to correspond to assigned BC.
        self.wink = pylith3d.allocateDouble(
            self.nwink)
        self.iwink = pylith3d.allocateInt(
            2*self.nwink)

        pylith3d.assign_wink(
            self.winkdef,
            self.wink,
            self.iwinkdef,
            self.iwinkid,
            self.iwink,
            self.id,
            self.numnp,
            self.nwink,
            self.nwinke)

        # Assign appropriate equation numbers to Iwinkx array, and compact Winkx
        # array to correspond to assigned BC.
        self.winkx = pylith3d.allocateDouble(
            self.nwinkx)
        self.iwinkx = pylith3d.allocateInt(
            2*self.nwinkx)

        pylith3d.assign_wink(
            self.winkxdef,
            self.winkx,
            self.iwinkxdef,
            self.iwinkxid,
            self.iwinkx,
            self.idx,
            self.numnp,
            self.nwinkx,
            self.nwinkxe)

        self.trace.log("Hello from pl3dsetup.numberequations (end)!")
            
        return


    def nfind(self):
        # Temporary arrays
        xtmp = pylith3d.allocateDouble(
            self.numsn)
        itmp = pylith3d.allocateInt(
            self.numsn)
        itmp1 = pylith3d.allocateInt(
            self.numsn)
        itmp2 = pylith3d.allocateInt(
            self.numsn)

        pylith3d.nfind(
            self.x,
            xtmp,
            self.idslp,
            self.ipslp,
            itmp,
            itmp1,
            itmp2,
            self.nslip,
            self.numslp,
            self.numsn,
            self.numnp)

        return


    def sortmesh(self):

        # This function sorts elements into families and sorts all other items that are
        # affected by this.

        self.trace.log("Hello from pl3dsetup.sortmesh (begin)!")
        
        print "Renumbering elements, split nodes, and slippery nodes:"

        self.sort_elements()

        self.stateSize = self.elementSizeInfo[0]
        self.state0Size = self.elementSizeInfo[1]
        self.propertySize = self.elementSizeInfo[2]

        # Sort split node entries.
        pylith3d.sort_split_nodes(
            self.nfault,
            self.indxiel,
            self.numfn,
            self.numelv)

        # Sort slippery node entries.
        pylith3d.sort_slip_nodes(
            self.nslip,
            self.indxiel,
            self.numslp,
            self.numelv)
            
        self.trace.log("Hello from pl3dsetup.sortmesh (end)!")

        return


    def sort_elements(self):
        # Sort elements into families.  The sorted elements are contained
        # in array Iens, and the index array for the new ordering is
        # Indxiel.  The index array for the original ordering is Ielindx.
        # The original element node array (Ien) and the associated
        # material type array (Mat) may be deallocated after sorting.
        
        self.iens = pylith3d.allocateInt(
            self.nen*self.numelv)
        self.ivfamily = pylith3d.allocateInt(
            6*self.nvfamilies)

        self.indxiel = pylith3d.allocateInt(
            self.numelv)

        self.ielindx = pylith3d.allocateInt(
            self.numelv)

        ivftmp = pylith3d.allocateInt(
            self.nvfamilies)

	self.elementSizeInfo = pylith3d.sort_elements(
            self.ien,
            self.mat,
            self.infmatmod,
            self.ivflist,
            self.ivfamily,
            self.iens,
            ivftmp,
            self.indxiel,
            self.ielindx,
            self.nen,
            self.ngauss,
            self.maxvfamilies,
            self.nvfamilies,
            self.prestressFlag,
            self.numelv,
            self.numnp)
        
        self.ien = None ### DEALLOC
        self.mat = None ### DEALLOC
        self.ivflist = None ### DEALLOC

        return

        
    def sparsesetup(self):

        # This function sets up sparse matrix and associated storage.

        self.trace.log("Hello from pl3dsetup.sparsesetup (begin)!")
        
        print "Setting up sparse matrix storage:"
        
        self.autoprestrStage, \
        self.elasticStage, \
        self.viscousStage, \
        self.iterateEvent = pylith3d.setupPETScLogging()

        # Arrays to map element equation numbers to global
        # Localize global equation numbers in element index arrays.
        self.lm = pylith3d.allocateInt(
            constants.numberDegreesFreedom*self.connectivitySize)
        self.lmx = pylith3d.allocateInt(
            constants.numberDegreesFreedom*self.connectivitySize)
        self.lmf = pylith3d.allocateInt(
            self.connectivitySize)

        pylith3d.local(
            self.id,
            self.numnp,
            self.iens,
            self.lm,
            self.numelv,
            self.nen)

        pylith3d.localf(
            self.iens,
            self.lmf,
            self.numelv,
            self.nfault,
            self.numfn,
            self.nen)

        pylith3d.localx(
            self.idx,
            self.numnp,
            self.iens,
            self.lmx,
            self.numelv,
            self.nslip,
            self.numslp,
            self.nen)

        # Keeping this for now as it may be wanted for output
        # self.nslip = None ### DEALLOC

        # Allocate and populate sparse matrix arrays.  Some of these are
        # temporary and are then deleted after use.
        iwork = pylith3d.cmp_stiffsz(
            self.neq,
            self.lm,
            self.lmx,
            self.numelv,
            self.numsn,
            self.nen)

        # Temporary arrays
        indx = pylith3d.allocateInt(
            self.neq)
        link = pylith3d.allocateInt(
            iwork)
        nbrs = pylith3d.allocateInt(
            iwork)

        stiffnessMatrixInfo = pylith3d.lnklst(
            self.neq,
            self.lm,
            self.lmx,
            self.numelv,
            self.nen,
            self.nee,
            indx,
            link,
            nbrs,
            iwork,
            self.numsn)

        self.nnz = stiffnessMatrixInfo[0]

        self.A, self.rhs, self.sol = pylith3d.createPETScMat(self.mesh)

        stiffnessMatrixStats = pylith3d.makemsr(
            self.A,
            indx,
            link,
            nbrs,
            self.neq,
            self.nnz,
            iwork)

        self.nmin = stiffnessMatrixStats[0]
        self.nmax = stiffnessMatrixStats[1]
        self.wavg = float(stiffnessMatrixStats[2])

	print ""
	print ""
        print "Sparse matrix information:"
	print ""
        print "numberGlobalEquations:     %i" % self.neq
        print "workingArraySize:          %i" % iwork
        print "stiffnessMatrixSize:       %i" % (self.nnz-1)
        print "stiffnessOffDiagonalSize:  %i" % stiffnessMatrixInfo[1]
        print "minimumNonzeroTermsPerRow: %i" % self.nmin
        print "maximumNonzeroTermsPerRow: %i" % self.nmax
        print "averageNonzeroTermsPerRow: %g" % self.wavg
	print ""
        
        self.trace.log("Hello from pl3dsetup.sparsesetup (end)!")

        return
        
    def allocateremaining(self):

        # This function allocates all remaining arrays that are needed for computations.
        
        self.trace.log("Hello from pl3dsetup.allocateremaining (begin)!")
        
        print "Allocating remaining storage:"
        
        # Create necessary lists and convert them to arrays
        self.grav = pylith3d.doubleListToArray(
            [self.gravityX.value,
             self.gravityY.value,
             self.gravityZ.value])
        
        # Allocate memory for all additional arrays

        # Force vectors
        if self.numtractions != 0:
            tractionFlag = 1
        else:
            tractionFlag = 0
        if self.gravityX.value != 0.0 or self.gravityY.value != 0.0 or self.gravityZ.value != 0.0:
            gravityFlag = 1
        else:
            gravityFlag = 0
        if self.numberConcForces != 0 or self.numdif != 0:
            concForceFlag = 1
        else:
            concForceFlag = 0
        if tractionFlag != 0 or gravityFlag != 0 or concForceFlag != 0:
            externFlag = 1
        else:
            externFlag = 0
	if self.nwink != 0:
	    winklerFlag = 1
        else:
            winklerFlag = 0
	if self.nwinkx != 0:
	    slipperyWinklerFlag = 1
        else:
            slipperyWinklerFlag = 0

        self.bextern = pylith3d.allocateDouble(
            externFlag*self.neq)
        self.btraction = pylith3d.allocateDouble(
            tractionFlag*self.neq)
        self.bgravity = pylith3d.allocateDouble(
            gravityFlag*self.neq)
        self.bconcForce = pylith3d.allocateDouble(
            concForceFlag*self.neq)
        self.bwink = pylith3d.allocateDouble(
            winklerFlag*self.neq)
        self.bwinkx = pylith3d.allocateDouble(
            slipperyWinklerFlag*self.neq)
        self.bintern = pylith3d.allocateDouble(
            self.neq)
        self.bresid = pylith3d.allocateDouble(
            self.neq)
        self.dispVec = pylith3d.allocateDouble(
            self.neq)
        self.dprev = pylith3d.allocateDouble(
            self.nprevdflag*self.neq)
            
        # Displacement arrays
        self.d = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numnp)
        self.deld = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numnp)
        self.dcur = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numnp)

        # Slippery node arrays
        self.dx = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numnp)
        self.deldx = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numnp)
        self.dxcur = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numnp)

        # Split node arrays
        self.dfault = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numfn)
        self.tfault = pylith3d.allocateDouble(
            constants.numberDegreesFreedom*self.numfn)

        # Local stiffness matrix arrays
        self.s = pylith3d.allocateDouble(
            self.maxElementEquations*self.maxElementEquations)
        self.stemp = pylith3d.allocateDouble(
            self.maxElementEquations*self.maxElementEquations)

        # Element arrays
        self.state = pylith3d.allocateDouble(
            self.stateSize)
        self.dstate = pylith3d.allocateDouble(
            self.stateSize)
        self.dmat = pylith3d.allocateDouble(
            constants.materialMatrixDimension*
            self.ngauss*
            self.numelv)
        self.iddmat = pylith3d.intListToArray( 
            constants.listIddmat)
        self.state0 = pylith3d.allocateDouble(
            self.state0Size)

        # Create arrays from lists that will be needed for the solution

        # nforce array
        self.nforce = pylith3d.intListToArray(
            [externFlag,
             tractionFlag,
             gravityFlag,
             concForceFlag,
             self.prestressFlag,
             winklerFlag,
             slipperyWinklerFlag,
             self.nprevdflag])
           
        # ncodat array
        self.ncodat = pylith3d.intListToArray(
            [self.icode,
             self.idebug])
            
        # npar array
        self.npar = pylith3d.intListToArray(
            [self.numelv,
             self.numat,
             self.numtractions,
             self.numslp,
             self.numfn,
             self.ipstrs,
             self.ipauto,
             self.stateSize,
             self.state0Size,
             self.nvfamilies,
             self.numdif,
             self.intord])

        # nprint array
        self.nprint = pylith3d.intListToArray(
            [self.icontr,
             self.idout,
             self.idsk,
             self.iucd])

        # nsysdat array
        self.nsysdat = pylith3d.intListToArray(
            [self.numnp,
             self.neq,
             self.nnz,
             self.numrot,
             self.numberPrestressEntries,
             self.numsn,
             self.numflt,
             self.propertySize,
             self.nwink,
             self.nwinkx,
             self.iskopt])

        # nunits array
        self.nunits = pylith3d.intListToArray(
            [self.f77StandardInput,
             self.f77StandardOutput,
             self.kr,
             self.kw,
             self.kp,
             self.kucd])

        # nvisdat array
        self.nvisdat = pylith3d.intListToArray(
            [self.ncycle,
             self.nintg,
             self.lastep,
             self.nhist])
        
        # rgiter array
        self.rgiter = pylith3d.doubleListToArray(
            [self.stol.value,
             self.dtol,
             self.epert])
        
        # rtimdat array
        self.rtimdat = pylith3d.doubleListToArray(
            [0.0, # currentTimeStepSize
             0.0, # currentAlfaParameter
             self.tpois,
             self.tyoungs.value])

        # ntimdat array
        self.ntimdat = pylith3d.intListToArray(
            [0, # currentTimeStep
             0, # currentIterationsBetweenReform
             0, # currentStepsBetweenReform
             0, # currentLargeDeformationFlag
             0, # currentMaximumIterations
             0, # currentNumberTotalIterations
             0, # currentNumberReforms
             0, # currentNumberTotalPcgIterations
             0, # reformFlagInt
             ])

        self.trace.log("Hello from pl3dsetup.allocateremaining (end)!")

        return


    def meshwrite(self):

        # This function outputs mesh information.
        # In the near future, this needs to be broken into classes for
        # Ascii output, plot output, UCD output, etc.

        self.trace.log("Hello from pl3dsetup.meshwriteascii (begin)!")
        
        print "Outputting Ascii mesh information:"

        # Write out global parameters
        pylith3d.write_global_info(
            self.title,
            self.idout,
            self.idsk,
            self.numnp,
            self.icode,
            self.idebug,
            self.kw,
            self.kp,
            self.ofile,
            self.pfile)

        # Write out nodal coordinates
        pylith3d.write_coords(
            self.x,
            self.numnp,
            self.kw,
            self.kp,
            self.idout,
            self.idsk,
            self.ofile,
            self.pfile)

        # Write out nodal boundary condition info
        pylith3d.write_bc(
            self.bond,
            self.ibond,
            self.numnp,
            self.kw,
            self.idout,
            self.ofile)

        # Write out local coordinate rotations
        pylith3d.write_skew(
            self.skew,
            self.numrot,
            self.iskopt,
            self.numnp,
            self.kw,
            self.idout,
            self.ofile)

        # Write stress computation and subiteration parameters.
        pylith3d.write_strscomp(
            self.stol.value,
            self.dtol,
            self.epert,
            self.kw,
            self.idout,
            self.ofile)

        pylith3d.write_subiter(
            self.nprevdflag,
            self.kw,
            self.idout,
            self.ofile)

        # Write out time step information
        pylith3d.write_timdat(
            self.delt,
            self.alfa,
            self.utol,
            self.ftol,
            self.etol,
            self.times,
            self.maxstp,
            self.maxit,
            self.ntdinit,
            self.lgdef,
            self.itmax,
            self.nintg,
            self.lastep,
            self.kw,
            self.idout,
            self.ofile)

        # Write out timesteps when full output is desired
        pylith3d.write_fuldat(
            self.iprint,
            self.icontr,
            self.icode,
            self.ncycle,
            self.lastep,
            self.kw,
            self.kp,
            self.idout,
            self.idsk,
            self.ofile,
            self.pfile)

        # Write out state variables desired for output
        pylith3d.write_stateout(
            self.istatout,
            self.nstatout,
            self.kw,
            self.kp,
            self.idout,
            self.idsk,
            self.ofile,
            self.pfile)

        # Write out load history information and deallocate Times array
        pylith3d.write_hist(
            self.histry,
            self.times,
            self.nhist,
            self.lastep,
            self.kw,
            self.idout,
            self.ofile)

        self.times = None ### DEALLOC

        # Write element info
        pylith3d.write_element_info(
            self.numelv,
            self.nen,
	    self.ngauss,
            self.ietypev,
            self.intord,
            self.ipstrs,
            self.ipauto,
            self.tpois,
            self.tyoungs.value,
            self.kw,
            self.idout,
            self.ofile)

        # Write element node array and deallocate Indxiel
        pylith3d.write_connect(
            self.iens,
            self.ivfamily,
            self.indxiel,
            self.nen,
	    self.ngauss,
            self.numelv,
            self.ietypev,
            self.nvfamilies,
            self.kw,
            self.kp,
            self.idout,
            self.idsk,
            self.ofile,
            self.pfile)

        self.indxiel = None ### DEALLOC

        # Write material properties
        pylith3d.write_props(
            self.prop,
            self.grav,
            self.ivfamily,
            self.infmatmod,
            self.nvfamilies,
            self.propertySize,
            self.idout,
            self.idsk,
            self.kw,
            self.kp,
            self.ofile,
            self.pfile)

        # Write mesh info to UCD file, if requested
        if self.iucd >= 0:
            pylith3d.write_ucd_mesh(
                self.x,
                self.numnp,
                self.iens,
                self.ivfamily,
                self.numelv,
                self.nvfamilies,
                self.sh,
                self.nen,
                self.ngauss,
                self.ietypev,
                self.istatout,
                self.nstatout,
                self.kucd,
                self.iucd,
                self.ucdroot)

        # Write traction info
        pylith3d.write_tractions(
            self.tractionverts,
            self.tractionvals,
            self.numtractions,
            self.nsnodes,
            self.kw,
            self.idout,
            self.ofile)
    
        # Write split node info
        pylith3d.write_split(
            self.fault,
            self.nfault,
            self.numfn,
            self.kw,
            self.kp,
            self.idout,
            self.idsk,
            self.ofile,
            self.pfile)

        # Write slippery node info
        pylith3d.write_slip(
            self.nslip,
            self.numslp,
            self.numsn,
            self.kw,
            self.kp,
            self.idout,
            self.idsk,
            self.ofile,
            self.pfile)

        # Write differential force info and deallocate Nslip
        pylith3d.write_diff(
            self.diforc,
            self.nslip,
            self.idhist,
            self.numslp,
            self.numdif,
            self.numnp,
            self.kw,
            self.idout,
            self.ofile)

        self.nslip = None ### DEALLOC

        # Write split nodes to plot file, if requested and deallocate Idftn
        pylith3d.write_split_plot(
            self.idftn,
            self.numflt,
            self.kp,
            self.idsk,
            self.pfile)

        self.idftn = None ### DEALLOC

        # Write Winkler force info and deallocate definition arrays
        pylith3d.write_wink(
            self.winkdef,
            self.iwinkdef,
            self.iwinkid,
            self.nwinke,
            self.kw,
            self.idout,
            self.ofile)

        self.winkdef = None ### DEALLOC
        self.iwinkdef = None ### DEALLOC

        # Write slippery node Winkler force info and deallocate definition arrays
        pylith3d.write_winkx(
            self.winkxdef,
            self.iwinkxdef,
            self.iwinkxid,
            self.nwinkxe,
            self.kw,
            self.idout,
            self.ofile)

        self.winkxdef = None ### DEALLOC
        self.iwinkxdef = None ### DEALLOC

        # Write sparse matrix info
        pylith3d.write_sparse_info(
            self.neq,
            self.nnz,
            self.nmin,
            self.nmax,
            self.wavg,
            self.idout,
            self.kw,
            self.ofile)

        self.trace.log("Hello from pl3dsetup.meshwrite (end)!")

        return


# The function of this code is to call the elastic and time-dependent solution
# drivers.  To do this, a number of previously-defined parameters need to be
# bundled into lists.


    def solveElastic(self):
        pylith3d.elastc(
            self.A,                # sparse
            self.rhs,
            self.sol,
            self.bextern,          # force
            self.btraction,
            self.bgravity,
            self.bconcForce,
            self.bintern,
            self.bresid,
            self.bwink,
            self.bwinkx,
            self.dispVec,
            self.dprev,
            self.nforce,
            self.grav,
            self.x,                # global
            self.d,
            self.deld,
            self.dcur,
            self.id,
            self.iwink,
            self.wink,
            self.nsysdat,
            self.iddmat,
            self.ibond,            # BC
            self.bond,
            self.dx,               # slip
            self.deldx,
            self.dxcur,
            self.diforc,
            self.idx,
            self.iwinkx,
            self.winkx,
            self.idslp,
            self.ipslp,
            self.idhist,
            self.fault,            # fault
            self.nfault,
            self.dfault,
            self.tfault,
            self.s,                # stiff
            self.stemp,
            self.state,            # element
            self.dstate,
            self.state0,
            self.dmat,
            self.iens,
            self.lm,
            self.lmx,
            self.lmf,
            self.ivfamily,
            self.npar,
            self.ielindx,
            self.tractionverts,    # traction
            self.tractionvals,
            self.gauss2d,
            self.sh2d,
            self.infetype2d,
            self.prop,             # material
            self.infmatmod,
            self.gauss,            # eltype
            self.sh,
            self.shj,
            self.infetype,
            self.histry,           # timdat
            self.rtimdat,
            self.ntimdat,
            self.nvisdat,
            self.maxstp,
            self.delt,
            self.alfa,
            self.maxit,
            self.ntdinit,
            self.lgdef,
            self.utol,
            self.ftol,
            self.etol,
            self.itmax,
            self.rgiter,           # stresscmp
            self.skew,             # skew
            self.ncodat,           # ioinfo
            self.nunits,
            self.nprint,
            self.istatout,
            self.nstatout,
            self.ofile,            # files
            self.pfile,
            self.ucdroot,
            self.elasticStage,     # PETSc logging
            self.iterateEvent)
        return

    def interpolatePoints(self, points):
        return pylith3d.interpolatePoints(self.mesh, self.sol, points)

    def runSimulation(self):
        # First define all of the lists that maintain variable values.  The
        # variables in these lists are altered during the running of the code
        # and should not be accessed directly except as a member of the list.
        # They should not have been defined previously.

        self.trace.log("Hello from pl3drun.run (begin)!")
        
        print "Beginning problem solution:"

        if False: # Temporarily out-of-order
            # Output approximate memory usage
            self.memorySizeMB =0.0
            self.memorySizeMB=self.memorySize/(1024.0*1024.0)

            print ""
            print "Approximate memory allocation for f77 arrays (MB): %g" % self.memorySizeMB

        # print "Just before pylith3d.autoprestr:"

        # Compute gravitational prestresses, if requested.
        if self.analysisType == "elasticSolution" or self.analysisType == "fullSolution":
            if self.ipstrs == 1:
                pylith3d.autoprestr(
                    self.A,                # sparse
                    self.rhs,
                    self.sol,
                    self.bextern,          # force
                    self.btraction,
                    self.bgravity,
                    self.bconcForce,
                    self.bintern,
                    self.bresid,
                    self.bwink,
                    self.bwinkx,
                    self.dispVec,
                    self.dprev,
                    self.nforce,
                    self.grav,
                    self.x,                # global
                    self.d,
                    self.deld,
                    self.dcur,
                    self.id,
                    self.iwink,
                    self.wink,
                    self.nsysdat,
                    self.iddmat,
                    self.ibond,            # BC
                    self.bond,
                    self.dx,               # slip
                    self.deldx,
                    self.dxcur,
                    self.diforc,
                    self.idx,
                    self.iwinkx,
                    self.winkx,
                    self.idslp,
                    self.ipslp,
                    self.idhist,
                    self.fault,            # split
                    self.nfault,
                    self.dfault,
                    self.tfault,
                    self.s,                # stiff
                    self.stemp,
                    self.state,            # element
                    self.dstate,
                    self.state0,
                    self.dmat,
                    self.iens,
                    self.lm,
                    self.lmx,
                    self.lmf,
                    self.ivfamily,
                    self.npar,
                    self.ielindx,
                    self.tractionverts,    # traction
                    self.tractionvals,
                    self.gauss2d,
                    self.sh2d,
                    self.infetype2d,
                    self.prop,             # material
                    self.infmatmod,
                    self.gauss,            # eltype
                    self.sh,
                    self.shj,
                    self.infetype,
                    self.histry,           # timdat
                    self.rtimdat,
                    self.ntimdat,
                    self.nvisdat,
                    self.maxstp,
                    self.delt,
                    self.alfa,
                    self.maxit,
                    self.ntdinit,
                    self.lgdef,
                    self.utol,
                    self.ftol,
                    self.etol,
                    self.itmax,
                    self.rgiter,           # stresscmp
                    self.skew,             # skew
                    self.ncodat,           # ioinfo
                    self.nunits,
                    self.nprint,
                    self.istatout,
                    self.nstatout,
                    self.ofile,            # files
                    self.pfile,
                    self.ucdroot,
                    self.autoprestrStage,  # PETSc logging
                    self.iterateEvent)

            # Perform elastic solution, if requested.
            self.solveElastic()
            pylith3d.outputMesh(self.fileRoot, self.mesh, self.sol)

        # Perform time-dependent solution, if requested.

        if self.analysisType == "fullSolution" and self.nintg > 1:
            if self.pythonTimestep:
                # Setup timestepping
                #   Open output files
                pylith3d.viscos_setup(self.nprint,
                                      self.nunits,
                                      self.ofile,
                                      self.pfile,
                                      self.viscousStage)
                numCycles         = pylith3d.intListRef(self.nvisdat, 0)
                numTimeStepGroups = pylith3d.intListRef(self.nvisdat, 1)
                numslp            = pylith3d.intListRef(self.npar, 3)
                iskopt            = pylith3d.intListRef(self.nsysdat, 10)
                icontr            = pylith3d.intListRef(self.nprint, 0)
                indexx            = 1 # Fortran index
                totalSteps        = 0 # This is ntot
                for cycle in range(numCycles):
                    if numCycles > 1: print '     working on cycle %d' % cycle
                    nextStartStep = 0 # This is naxstp
                    timeStep      = 0 # This is nstep
                    startStep     = 0 # This is nfirst
                    time          = 0.0

                    for tsGroup in range(1, numTimeStepGroups):
                        # Define constants
                        dt = pylith3d.doubleListRef(self.delt, tsGroup) # This is deltp
                        pylith3d.doubleListSet(self.rtimdat, 0, dt)
                        alfap = pylith3d.doubleListRef(self.alfa, tsGroup)
                        pylith3d.doubleListSet(self.rtimdat, 1, alfap)
                        pylith3d.intListSet(self.ntimdat, 0, timeStep)
                        maxitp = pylith3d.intListRef(self.maxit, tsGroup)
                        pylith3d.intListSet(self.ntimdat, 1, maxitp)
                        ntdinitp = pylith3d.intListRef(self.ntdinit, tsGroup)
                        pylith3d.intListSet(self.ntimdat, 2, ntdinitp)
                        lgdefp = pylith3d.intListRef(self.lgdef, tsGroup)
                        pylith3d.intListSet(self.ntimdat, 3, lgdefp)
                        itmaxp = pylith3d.intListRef(self.itmax, tsGroup)
                        pylith3d.intListSet(self.ntimdat, 4, itmaxp)
                        gtol = [pylith3d.doubleListRef(self.utol, tsGroup),
                                pylith3d.doubleListRef(self.ftol, tsGroup),
                                pylith3d.doubleListRef(self.etol, tsGroup)]
                        startStep     = nextStartStep + 1
                        nextStartStep = startStep + pylith3d.intListRef(self.maxstp, tsGroup) - 1

                        ltim = 1

                        for j in range(startStep, nextStartStep+1):
                            totalSteps += 1
                            timeStep   += 1
                            pylith3d.intListSet(self.ntimdat, 0, timeStep)
                            time += dt
                            skc   = (numslp != 0 and (iskopt == 2 or (iskopt <= 0 and abs(iskopt) == timeStep)))

                            pylith3d.viscos_step(
                                self.A,                # sparse
                                self.rhs,
                                self.sol,
                                self.bextern,          # force
                                self.btraction,
                                self.bgravity,
                                self.bconcForce,
                                self.bintern,
                                self.bresid,
                                self.bwink,
                                self.bwinkx,
                                self.dispVec,
                                self.dprev,
                                self.nforce,
                                self.grav,
                                self.x,                # global
                                self.d,
                                self.deld,
                                self.dcur,
                                self.id,
                                self.iwink,
                                self.wink,
                                self.nsysdat,
                                self.iddmat,
                                self.ibond,            # BC
                                self.bond,
                                self.dx,               # slip
                                self.deldx,
                                self.dxcur,
                                self.diforc,
                                self.idx,
                                self.iwinkx,
                                self.winkx,
                                self.idslp,
                                self.ipslp,
                                self.idhist,
                                self.fault,            # fault
                                self.nfault,
                                self.dfault,
                                self.tfault,
                                self.s,                # stiff
                                self.stemp,
                                self.state,            # element
                                self.dstate,
                                self.state0,
                                self.dmat,
                                self.iens,
                                self.lm,
                                self.lmx,
                                self.lmf,
                                self.ivfamily,
                                self.npar,
                                self.ielindx,
                                self.tractionverts,    # traction
                                self.tractionvals,
                                self.gauss2d,
                                self.sh2d,
                                self.infetype2d,
                                self.prop,             # material
                                self.infmatmod,
                                self.gauss,            # eltype
                                self.sh,
                                self.shj,
                                self.infetype,
                                self.histry,           # timdat
                                self.rtimdat,
                                self.ntimdat,
                                self.nvisdat,
                                self.maxstp,
                                self.delt,
                                self.alfa,
                                self.maxit,
                                self.ntdinit,
                                self.lgdef,
                                self.utol,
                                self.ftol,
                                self.etol,
                                self.itmax,
                                self.rgiter,           # stresscmp
                                self.skew,             # skew
                                self.iprint,           # ioinfo
                                self.ncodat,
                                self.nunits,
                                self.nprint,
                                self.istatout,
                                self.nstatout,
                                self.ofile,            # files
                                self.pfile,
                                self.ucdroot,
                                self.viscousStage,     # PETSc logging
                                self.iterateEvent,
                                totalSteps,
                                ltim,
                                indexx,
                                cycle,
                                tsGroup,
                                j,
                                skc,
                                startStep,
                                timeStep,
                                time,
                                dt,
                                lgdefp,
                                gtol)
                            ltim = 0
                            if (totalSteps == pylith3d.intListRef(self.iprint, indexx-1)):
                                pylith3d.outputMesh(self.fileRoot+'.'+str(totalSteps), self.mesh, self.sol)
                                indexx += 1
                            if (indexx > icontr): indexx = icontr
                print " Total number of equilibrium iterations        =",pylith3d.intListRef(self.ntimdat, 5)
                print " Total number of stiffness matrix reformations =",pylith3d.intListRef(self.ntimdat, 6)
                print " Total number of displacement subiterations    =",pylith3d.intListRef(self.ntimdat, 7)
                pylith3d.viscos_cleanup(self.ntimdat, self.nprint, self.nunits)
            else:
                pylith3d.viscos(
                    self.A,                # sparse
                    self.rhs,
                    self.sol,
                    self.bextern,          # force
                    self.btraction,
                    self.bgravity,
                    self.bconcForce,
                    self.bintern,
                    self.bresid,
                    self.bwink,
                    self.bwinkx,
                    self.dispVec,
                    self.dprev,
                    self.nforce,
                    self.grav,
                    self.x,                # global
                    self.d,
                    self.deld,
                    self.dcur,
                    self.id,
                    self.iwink,
                    self.wink,
                    self.nsysdat,
                    self.iddmat,
                    self.ibond,            # BC
                    self.bond,
                    self.dx,               # slip
                    self.deldx,
                    self.dxcur,
                    self.diforc,
                    self.idx,
                    self.iwinkx,
                    self.winkx,
                    self.idslp,
                    self.ipslp,
                    self.idhist,
                    self.fault,            # fault
                    self.nfault,
                    self.dfault,
                    self.tfault,
                    self.s,                # stiff
                    self.stemp,
                    self.state,            # element
                    self.dstate,
                    self.state0,
                    self.dmat,
                    self.iens,
                    self.lm,
                    self.lmx,
                    self.lmf,
                    self.ivfamily,
                    self.npar,
                    self.ielindx,
                    self.tractionverts,    # traction
                    self.tractionvals,
                    self.gauss2d,
                    self.sh2d,
                    self.infetype2d,
                    self.prop,             # material
                    self.infmatmod,
                    self.gauss,            # eltype
                    self.sh,
                    self.shj,
                    self.infetype,
                    self.histry,           # timdat
                    self.rtimdat,
                    self.ntimdat,
                    self.nvisdat,
                    self.maxstp,
                    self.delt,
                    self.alfa,
                    self.maxit,
                    self.ntdinit,
                    self.lgdef,
                    self.utol,
                    self.ftol,
                    self.etol,
                    self.itmax,
                    self.rgiter,           # stresscmp
                    self.skew,             # skew
                    self.iprint,           # ioinfo
                    self.ncodat,
                    self.nunits,
                    self.nprint,
                    self.istatout,
                    self.nstatout,
                    self.ofile,            # files
                    self.pfile,
                    self.ucdroot,
                    self.viscousStage,     # PETSc logging
                    self.iterateEvent)
        pylith3d.destroyPETScMat(self.A,self.rhs,self.sol)

        self.trace.log("Hello from pl3drun.run (end)!")
        
        return


# end of file 
