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


# simple C external declarations
cimport constants
cimport elementtypes
cimport libpylith3d
cimport petsc

# external Pyrex modules
cimport PyLithMeshLib

# Pyrex code inlined in this module
include "array.pyx"
include "exceptionhandler.pxd"
include "f77io.pyx"
include "setup.pyx"


cdef enum:
    prestress = 0 # code for reading prestress input files is presently disabled


cdef class PyLith:

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # Public Interface
    #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # These properties and methods are accessible from Python.  They
    # form the interface between the UI (the Pyre Component) and the
    # native code.

    cdef public object title

    # Basename for all files (may be overridden by specific filename entries).
    cdef public object fileRoot

    # Output filenames (all are optional).
    cdef public object ofile      # asciiOutputFile    ooo.ascii
    cdef public object pfile      # plotOutputFile     ooo.plot
    cdef public object ucdroot    # ucdOutputRoot      ooo.*.inp

    # Required input files.
    cdef public object coordinateInputFile            # coordinateInputFile            iii.coord
    cdef public object bcfile                         # bcInputFile                    iii.bc
    cdef public object timeStepInputFile              # timeStepInputFile              iii.time
    cdef public object stfile                         # stateVariableInputFile         iii.statevar
    cdef public object connectivityInputFile          # connectivityInputFile          iii.connect

    # This file is only required for time-dependent problems.
    cdef public object fofile                         # fullOutputInputFile            iii.fuldat

    # These files are optional unless generating Green's functions, in
    # which case they are required.
    cdef public object sampleLocationFile             # sampleLocationFile             iii.sample
    cdef public object spfile                         # splitNodeInputFile             iii.split

    # Optional input files.
    cdef public object skfile                         # rotationInputFile              iii.skew
    cdef public object hfile                          # loadHistoryInputFile           iii.hist
    cdef public object tractionInputFile              # tractionInputFile              iii.traction

    # Unused input files.
    cdef public object wfile                          # winklerInputFile               iii.wink
    cdef public object materialHistoryInputFile       # materialHistoryInputFile       iii.mhist
    cdef public object prestressInputFile             # prestressInputFile             iii.prestr
    cdef public object slfile                         # slipperyNodeInputFile          iii.slip
    cdef public object difile                         # differentialForceInputFile     iii.diff
    cdef public object wxfile                         # slipperyWinklerInputFile       iii.winkx

    # Output option flags.
    cdef public int idout             # asciiOutput
    cdef public int idsk              # plotOutput
    cdef public int iucd              # ucdOutput

    # Additional option flags.
    cdef public int icode             # analysisType
    cdef public int pythonTimestep    # pythonTimestep
    cdef public int generateGreen     # generateGreen
    cdef public int idebug            # debuggingOutput
    cdef public int ncycle            # numberCycles

    # Unused option flags.
    cdef public int iskopt            # autoRotateSlipperyNodes

    #
    # Category 2 parameters formerly placed in *.keyval files.
    #

    cdef public double stol           # stressTolerance in Pa
    cdef public double dtol           # minimumStrainPerturbation
    cdef public double epert          # initialStrainPerturbation

    cdef public int    nprevdflag     # usePreviousDisplacementFlag

    cdef public int    intord         # quadratureOrder

    cdef public int    ipstrs         # prestressAutoCompute
    cdef public int    ipauto         # prestressAutoChangeElasticProps
    cdef public double tpois          # prestressAutoComputePoisson
    cdef public double tyoungs        # prestressAutoComputeYoungs in Pa

    #
    # Array properties.
    #

    property grav: # gravityX, gravityY, gravityZ in m/(s*s)
        def __set__(self, gravity):
            cdef int dim
            dim = sizeof(self._grav)/sizeof(self._grav[0])
            assert(len(gravity) == dim)
            for i from 0 <= i < dim:
                self._grav[i] = gravity[i]

    property prscal: # prestressScale{Xx,Yy,Zz,Xy,Xz,Yz}
        def __set__(self, prestressScale):
            cdef int dim
            dim = sizeof(self._prscal)/sizeof(self._prscal[0])
            assert(len(prestressScale) == dim)
            for i from 0 <= i < dim:
                self._prscal[i] = prestressScale[i]

    property wscal: # winklerScaleX, winklerScaleY, winklerScaleZ
        def __set__(self, winklerScale):
            cdef int dim
            dim = sizeof(self._wscal)/sizeof(self._wscal[0])
            assert(len(winklerScale) == dim)
            for i from 0 <= i < dim:
                self._wscal[i] = winklerScale[i]
    
    property wxscal: # winklerSlipScaleX, winklerSlipScaleY, winklerSlipScaleZ
        def __set__(self, winklerSlipScale):
            cdef int dim
            dim = sizeof(self._wxscal)/sizeof(self._wxscal[0])
            assert(len(winklerSlipScale) == dim)
            for i from 0 <= i < dim:
                self._wxscal[i] = winklerSlipScale[i]

    #
    # materialPropertiesInputFile    iii.prop
    #

    cdef public int   numat
    
    property prop:
        def __set__(self, propertyList):
            cdef int nprops
            nprops = len(propertyList)
            self._prop = DoubleArray(nprops)
            for i from 0 <= i < nprops:
                self._prop.ptr[i] = propertyList[i]

    property infmat:
        def __set__(self, materialModel):
            cdef int ninfmats
            ninfmats = len(materialModel)
            self._infmat = IntArray(ninfmats)
            for i from 0 <= i < ninfmats:
                self._infmat.ptr[i] = materialModel[i]


    # Unit numbers for Fortran I/O.
    cdef public int   kti             # f77StandardInput
    cdef public int   kto             # f77StandardOutput
    cdef public int   kr              # f77FileInput
    cdef public int   kw              # f77AsciiOutput
    cdef public int   kp              # f77PlotOutput
    cdef public int   kucd            # f77UcdOutput

    #
    # Mesh.
    #

    ## work-around silly Pyrex limitation
    #cdef public PyLithMeshLib.Mesh mesh
    cdef PyLithMeshLib.Mesh _mesh
    property mesh:
        def __get__(self): return self._mesh
        def __set__(self, mesh): self._mesh = mesh


    def run(self, points):

        self.scan()
        self.read()
        self.numberequations()
        self.sortmesh()
        self.sparsesetup()
        self.allocateremaining()

        #
        # open output files
        #

        # xxx.ascii
        f77open(self.kw, self.ofile)

        # xxx.plot
        if self.idsk == 0:
            pass
        elif self.idsk == 1:
            f77open(self.kp, self.pfile)
        elif self.idsk == 2:
            f77open(self.kp, self.pfile, form="unformatted")
        else:
            raise ValueError("idsk (%d) is not in [0, 1, 2]" % idsk)

        try:
            self.meshwrite()

            if self.generateGreen:
                self.greenFunction(points)
            else:
                self.runSimulation()

        finally:
            # close output files
            f77close(self.kw)
            if self.idsk != 0: f77close(self.kp)

        return


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # Private
    #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # (public)
    cdef double _grav[3]
    cdef double _prscal[6]
    cdef double _wscal[3]
    cdef double _wxscal[3]
    cdef DoubleArray _prop
    cdef IntArray _infmat

    # debugging
    cdef object trace
    
    # (constants)
    cdef int maxElementEquations
    cdef int neni[10]
    
    # matmod_def
    cdef int infmatmod[6 * constants.maxMaterialModels]
    
    # scan_coords
    cdef int numnp
    cdef double cscale
    
    # scan_bc
    cdef int numbc
    cdef double dscale, vscale, fscale
    
    # scan_wink
    cdef int nwinke, nwink
    
    # scan_skew
    cdef int numrot
    cdef double runits
    
    # scan_timdat
    cdef int lastep, nintg
    cdef double tunits
    
    # scan_fuldat
    cdef int icontr

    # scan_hist
    cdef int nhist
    
    # scan_prestr
    cdef int nprestr
    cdef int nprestrflag
    
    # scan_tractions
    cdef int numtractions
    cdef double tscale
    
    # scan_split
    cdef int numfn
    
    # scan_slip
    cdef int numslp

    # scan_diff
    cdef int numdif
    
    # scan_winkx
    cdef int nwinkxe, nwinkx

    # scan_connect
    cdef int maxvfamilies
    cdef IntArray ivflist
    cdef int numelv, nvfamilies, ietypev

    # getdef, preshape, preshape2d
    cdef int nen, nsnodes, ngauss, nee
    cdef int connectivitySize
    cdef DoubleArray sh, shj, gauss
    cdef DoubleArray sh2d, gauss2d
    cdef int infetype[4], infetype2d[4]

    # read_coords
    cdef DoubleArray x

    # read_bc
    cdef DoubleArray bond
    cdef IntArray ibond
    cdef int numberConcForces

    # read_skew
    cdef DoubleArray skew

    # read_timdat
    cdef DoubleArray delt, alfa, utol, ftol, etol, times
    cdef IntArray maxstp, maxit, ntdinit, lgdef, itmax

    # read_fuldat
    cdef IntArray iprint

    # read_stateout
    cdef int istatout[3*constants.maxStateVariables]
    cdef int nstatout[3]

    # read_hist
    cdef DoubleArray histry

    # read_connect
    cdef IntArray ien, mat

    # read_tractions
    cdef IntArray tractionverts
    cdef DoubleArray tractionvals

    # read_split
    cdef IntArray nfault
    cdef DoubleArray fault
    cdef int numflt

    # read_slip
    cdef IntArray nslip
    cdef int numsn

    # read_diff
    cdef IntArray idhist
    cdef DoubleArray diforc

    # read_wink
    cdef IntArray iwinkdef, iwinkid
    cdef DoubleArray winkdef
    cdef IntArray iwinkxdef, iwinkxid
    cdef DoubleArray winkxdef

    # id_split
    cdef IntArray idftn

    # create_id
    cdef IntArray id, idx, idslp
    cdef int neq

    # nfind
    cdef IntArray ipslp

    # assign_wink
    cdef DoubleArray wink, winkx
    cdef IntArray iwink, iwinkx

    # sort_elements
    cdef IntArray iens, ivfamily, indxiel, ielindx
    cdef int nstatesz, nstatesz0, npropsz

    # setupPETScLogging
    cdef petsc.PetscInt autoprestrStage, elasticStage, viscousStage
    cdef petsc.PetscEvent iterateEvent

    # local, localf, localx
    cdef IntArray lm, lmf, lmx

    # lnklst
    cdef int nnz

    # makemsr
    cdef int nmin, nmax
    cdef double wavg

    # Force vectors
    cdef DoubleArray bextern, btraction, bgravity, bconcForce
    cdef DoubleArray bwink, bwinkx, bintern, bresid, dispVec, dprev
    
    # Displacement arrays
    cdef DoubleArray d, deld, dcur
    
    # Slippery node arrays
    cdef DoubleArray dx, deldx, dxcur

    # Split node arrays
    cdef DoubleArray dfault, tfault

    # Local stiffness matrix arrays
    cdef DoubleArray s, stemp

    # Element arrays
    cdef DoubleArray state, dstate, dmat, state0
    cdef int iddmat[36]

    # other arrays needed for the solution
    cdef int nforce[8]
    cdef int ncodat[2]
    cdef int npar[12]
    cdef int nprint[4]
    cdef int nsysdat[11]
    cdef int nunits[6]
    cdef int nvisdat[4]
    cdef double rgiter[3]
    cdef double rtimdat[4]
    cdef int ntimdat[9]
    

    def __new__(self):

        # Initialize public properties, just to be safe.  (These
        # defaults will be overwritten by the Pyre component.)
        
        # Output option flags.
        self.idout             = 1          # asciiOutput=echo
        self.idsk              = 0          # plotOutput=none
        self.iucd              = 0          # ucdOutput=none

        # Additional option flags.
        self.icode             = 3          # analysisType=fullSolution
        self.pythonTimestep    = 0          # pythonTimestep=False
        self.generateGreen     = 0          # generateGreen=False
        self.idebug            = 0          # debuggingOutput=False
        self.ncycle            = 1          # numberCycles=1

        # Unused option flags.
        self.iskopt            = 2          # autoRotateSlipperyNodes=True

        #
        # Category 2 parameters formerly placed in *.keyval files.
        #

        self.stol              = 1.0e-12    # stressTolerance
        self.dtol              = 1.0e-7     # minimumStrainPerturbation
        self.epert             = 1.0e-1     # initialStrainPerturbation

        self.nprevdflag        = 0          # usePreviousDisplacementFlag

        self.intord            = 1          # quadratureOrder=Full

        self.ipstrs            = 0          # prestressAutoCompute=False
        self.ipauto            = 0          # prestressAutoChangeElasticProps=False
        self.tpois             = 0.49       # prestressAutoComputePoisson
        self.tyoungs           = 1.0e30     # prestressAutoComputeYoungs in Pa

        # Array properties.

        cdef int dim
        
        # _grav
        dim = sizeof(self._grav)/sizeof(self._grav[0])
        for i from 0 <= i < dim:
            self._grav[i] = 0.0

        # _prscal
        dim = sizeof(self._prscal)/sizeof(self._prscal[0])
        for i from 0 <= i < dim:
            self._prscal[i] = 1.0

        # _wscal
        dim = sizeof(self._wscal)/sizeof(self._wscal[0])
        for i from 0 <= i < dim:
            self._wscal[i] = 1.0

        dim = sizeof(self._wxscal)/sizeof(self._wxscal[0])
        for i from 0 <= i < dim:
            self._wxscal[i] = 1.0

        # Unit numbers for Fortran I/O.
        self.kti               = 5          # f77StandardInput
        self.kto               = 6          # f77StandardOutput
        self.kr                = 10         # f77FileInput
        self.kw                = 11         # f77AsciiOutput
        self.kp                = 12         # f77PlotOutput
        self.kucd              = 13         # f77UcdOutput

        return


    def __init__(self):
        import journal
        self.trace = journal.debug("pylith3d.trace")


    cdef outputSampleValues(self, filename, impulseNode, values):
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


    cdef greenFunction(self, points):
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


# derived or automatically-specified quantities (category 3)

    cdef scan(self):

        self.trace.log("Hello from PyLith.scan (begin)!")
        self.trace.log("Scanning ascii files to determine dimensions.")

        import pyre.units
        
        cdef int errorcode
        cdef char errorstring[errorstring_size]
        errorcode = 0
        
        uparser = pyre.units.parser()

        # This is a test version where the geometry type is automatically
        # specified by using Pylith3d.  The geometry type is only used for
        # f77 routines and not in pyre. An integer value is also defined
        # for use in f77 routines.
        # Define some integer values that are derived from string variables.

        # Invariant parameters related to element type
        self.maxElementEquations = constants.numberDegreesFreedom*constants.maxElementNodes
        
        neni = [8, 7, 6, 5, 4, 20, 18, 15, 13, 10]
        cdef int dim
        dim = sizeof(self.neni)/sizeof(self.neni[0])
        assert(len(neni) == dim)
        for i from 0 <= i < dim:
            self.neni[i] = neni[i]

        # Invariant parameters related to material model
        libpylith3d.matmod_def(
            self.infmatmod  # intent(out)
            )

        # Parameters derived from the number of entries in a file.

        cdef char coord_units[30]
        libpylith3d.scan_coords(
            &self.numnp,  # intent(out)
            &self.kr,
            coord_units,  # intent(out)
            self.coordinateInputFile,
            &errorcode,
            errorstring,
            sizeof(coord_units),
            len(self.coordinateInputFile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)
        coord_units[sizeof(coord_units)-1] = 0   # null-terminate

        self.cscale = uparser.parse(str(coord_units).strip()).value

        cdef char displacement_units[30]
        cdef char velocity_units[30]
        cdef char force_units[30]
        libpylith3d.scan_bc(
            &self.numbc,          # intent(out)
            &self.kr,
            displacement_units,   # intent(out)
            velocity_units,       # intent(out)
            force_units,          # intent(out)
            self.bcfile,
            &errorcode,
            errorstring,
            sizeof(displacement_units),
            sizeof(velocity_units),
            sizeof(force_units),
            len(self.bcfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)
        displacement_units[sizeof(displacement_units)-1] = 0  # null-terminate
        velocity_units[sizeof(velocity_units)-1] = 0  # null-terminate
        force_units[sizeof(force_units)-1] = 0  # null-terminate
        
        if self.numbc > 0:
            self.dscale = uparser.parse(str(displacement_units).strip()).value
            self.vscale = uparser.parse(str(velocity_units).strip()).value
            self.fscale = uparser.parse(str(force_units).strip()).value
        else:
            self.dscale = 0.0
            self.vscale = 0.0
            self.fscale = 0.0

        libpylith3d.scan_wink(
            &self.nwinke,  # intent(out)
            &self.nwink,   # intent(out)
            &self.kr,
            self.wfile,
            &errorcode,
            errorstring,
            len(self.wfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        cdef char rotation_units[30]
        libpylith3d.scan_skew(
            &self.numrot,    # intent(out)
            &self.kr,
            rotation_units,  # intent(out)
            self.skfile,
            &errorcode,
            errorstring,
            sizeof(rotation_units),
            len(self.skfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)
        rotation_units[sizeof(rotation_units)-1] = 0  # null-terminate

        if self.numrot != 0:
            self.runits = uparser.parse(str(rotation_units).strip()).value
        else:
            self.runits = 0.0

        cdef char time_units[30]
        libpylith3d.scan_timdat(
            &self.lastep,  # intent(out)
            &self.nintg,   # intent(out)
            &self.kr,
            time_units,    # intent(out)
            self.timeStepInputFile,
            &errorcode,
            errorstring,
            sizeof(time_units),
            len(self.timeStepInputFile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)
        time_units[sizeof(time_units)-1] = 0  # null-terminate

        self.tunits = uparser.parse(str(time_units).strip()).value

        libpylith3d.scan_fuldat(
            &self.icode,
            &self.lastep,
            &self.icontr,  # intent(out)
            &self.kr,
            self.fofile,
            &errorcode,
            errorstring,
            len(self.fofile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        libpylith3d.scan_hist(
            &self.nhist,  # intent(out)
            &self.kr,
            self.hfile,
            &errorcode,
            errorstring,
            len(self.hfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.scan_connect()

        if prestress:
            ###### TODO
#             libpylith3d.scan_prestr(
#                 &constants.stateVariableDimension,
#                 &self.numberPrestressGaussPoints,
#                 &self.nprestr,  # intent(out)
#                 &self.numberElements,
#                 &self.ipstrs,
#                 &self.kr,
#                 self.prestressInputFile,
#                 &errorcode,
#                 errorstring,
#                 len(self.prestressInputFile),
#                 sizeof(errorstring))
            exceptionhandler(errorcode, errorstring)
            self.nprestr = 0
        else:
            self.nprestr = 0

        cdef char traction_units[30]
        cdef int nsnodesmax
        nsnodesmax = constants.nsnodesmax
        libpylith3d.scan_tractions(
            &self.numtractions,  # intent(out)
            &nsnodesmax,
            &self.kr,
            traction_units,      # intent(out)
            self.tractionInputFile,
            &errorcode,
            errorstring,
            sizeof(traction_units),
            len(self.tractionInputFile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)
        traction_units[sizeof(traction_units)-1] = 0  # null-terminate

        if self.numtractions != 0:
            self.tscale = uparser.parse(str(traction_units).strip()).value
        else:
            self.tscale = 0.0

        libpylith3d.scan_split(
            &self.numfn,  # intent(out)
            &self.kr,
            self.spfile,
            &errorcode,
            errorstring,
            len(self.spfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        libpylith3d.scan_slip(
            &self.numslp,  # intent(out)
            &self.kr,
            self.slfile,
            &errorcode,
            errorstring,
            len(self.slfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        libpylith3d.scan_diff(
            &self.numslp,
            &self.numdif,  # intent(out)
            &self.kr,
            self.difile,
            &errorcode,
            errorstring,
            len(self.difile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        libpylith3d.scan_winkx(
            &self.numslp,
            &self.nwinkxe,  # intent(out)
            &self.nwinkx,   # intent(out)
            &self.kr,
            self.wxfile,
            &errorcode,
            errorstring,
            len(self.wxfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.trace.log("Hello from PyLith.scan (end)!")

        return


    cdef scan_connect(self):
        
        cdef int errorcode
        cdef char errorstring[errorstring_size]
        errorcode = 0
        
        # At present, we assume that the number of element families is equal to
        # the number of material types used, since only one volume element type at a
        # time is allowed.
        self.maxvfamilies = (constants.numberAllowedVolumeElementTypes *
                             self.numat)

        self.ivflist = IntArray(3*self.maxvfamilies)

        libpylith3d.scan_connect(
            self.neni,
            self.infmatmod,
            self._infmat.ptr,
            self.ivflist.ptr,    # intent(out)
            &self.maxvfamilies,
	    &self.numat,
            &self.numelv,        # intent(out)
            &self.nvfamilies,    # intent(out)
            &self.ietypev,       # intent(out)
            &self.kr,
            self.connectivityInputFile,
            &errorcode,
            errorstring,
            len(self.connectivityInputFile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        return


# The main function of this code is to emulate the original functionality of
# input.f in the original version of TECTON.  This code segment controls the
# allocation of memory and the reading of the input file.  Additional functions
# covered by this code include the sparse matrix setup portion, which also does
# some memory allocation.  Additional code sections will call the main elastic
# and time-dependent solution drivers, which are presently f77 subroutines.


    cdef read(self):

        cdef int errorcode
        cdef char errorstring[errorstring_size]
        errorcode = 0
        
        # This function reads all input and performs some memory allocation.

        self.trace.log("Hello from pl3dsetup.read (begin)!")
        
        print "Reading problem definition and allocating necessary storage."

        # Set up global integration info.
        self.getdef()
        
        #
        # Node-based info (coordinates, displacement arrays, BC, and skew BC).
        #

        self.x = DoubleArray(constants.numberSpaceDimensions*self.numnp)
        libpylith3d.read_coords(
            self.x.ptr,      # intent(out)
            &self.cscale,
            &self.numnp,
            &self.kr,
            self.coordinateInputFile,
            &errorcode,
            errorstring,
            len(self.coordinateInputFile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.ibond = IntArray(constants.numberDegreesFreedom*self.numnp)
        self.bond = DoubleArray(constants.numberDegreesFreedom*self.numnp)
        libpylith3d.read_bc(
            self.bond.ptr,     # intent(out)
            &self.dscale,
            &self.vscale,
            &self.fscale,
            self.ibond.ptr,    # intent(out)
            &self.numnp,
            &self.numbc,
            &self.numberConcForces,    # intent(out)
            &self.kr,
            self.bcfile,
            &errorcode,
            errorstring,
            len(self.bcfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.skew = DoubleArray(constants.numberSkewDimensions*self.numnp)
        libpylith3d.read_skew(
            self.skew.ptr,    # intent(out)
            &self.runits,
            &self.numrot,
            &self.numnp,
            &self.iskopt,
            &self.kr,
            self.skfile,
            &errorcode,
            errorstring,
            len(self.skfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        #
        # Allocate and read time step, time output, and load history info.
        #

        self.delt = DoubleArray(self.nintg)
        self.alfa = DoubleArray(self.nintg)
        self.utol = DoubleArray(self.nintg)
        self.ftol = DoubleArray(self.nintg)
        self.etol = DoubleArray(self.nintg)
        # Note that array 'times' is needed for output, if requested.
        self.times = DoubleArray(self.lastep+1)
        self.maxstp = IntArray(self.nintg)
        self.maxit = IntArray(self.nintg)
        self.ntdinit = IntArray(self.nintg)
        self.lgdef = IntArray(self.nintg)
        self.itmax = IntArray(self.nintg)
        libpylith3d.read_timdat(
            self.delt.ptr,       # intent(out)
            self.alfa.ptr,       # intent(out)
            self.utol.ptr,       # intent(out)
            self.ftol.ptr,       # intent(out)
            self.etol.ptr,       # intent(out)
            self.times.ptr,      # intent(out)
            &self.tunits,
            self.maxstp.ptr,     # intent(out)
            self.maxit.ptr,      # intent(out)
            self.ntdinit.ptr,    # intent(out)
            self.lgdef.ptr,      # intent(out)
            self.itmax.ptr,      # intent(out)
            &self.nintg,
            &self.lastep,
            &self.kr,
            self.timeStepInputFile,
            &errorcode,
            errorstring,
            len(self.timeStepInputFile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.iprint = IntArray(self.icontr)
        libpylith3d.read_fuldat(
            self.iprint.ptr,    # intent(out)
            &self.icontr,
            &self.icode,
            &self.ncycle,
            &self.lastep,
            &self.kr,
            self.fofile,
            &errorcode,
            errorstring,
            len(self.fofile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        libpylith3d.read_stateout(
            self.istatout,
            self.nstatout,
            &self.kr,
            self.stfile,
            &errorcode,
            errorstring,
            len(self.stfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.histry = DoubleArray((self.lastep+1)*self.nhist)
        libpylith3d.read_hist(
            self.histry.ptr,    # intent(out)
            self.times.ptr,
            &self.nhist,
            &self.lastep,
            &self.kr,
            self.hfile,
            &errorcode,
            errorstring,
            len(self.hfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        #
        # Allocate and read info on connectivities and prestresses
        #
        
        if self.nprestr != 0 or self.ipstrs != 0:
            self.nprestrflag = 1
        else:
            self.nprestrflag = 0

        self.ien = IntArray(self.nen*self.numelv)
        self.mat = IntArray(self.numelv)
        libpylith3d.read_connect(
            self.ien.ptr,    # intent(out)
            self.mat.ptr,    # intent(out)
            &self.nen,
            &self.numelv,
            &self.numnp,
            &self.nvfamilies,
            &self.kr,
            self.connectivityInputFile,
            &errorcode,
            errorstring,
            len(self.connectivityInputFile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        if prestress:
            ###### TODO
#             libpylith3d.read_prestr(
#                 self.stn,
#                 self.st0,
#                 self._prscal,
#                 self.numberStressComponents,
#                 self.numberGaussPoints,
#                 self.numberPrestressGaussPoints,
#                 self.numberElements,
#                 self.nprestr,
#                 self.ipstrs,
#                 self.idout,
#                 self.kr,
#                 self.kw,
#                 self.prestressInputFile)
            pass

        # Read traction BC
        self.tractionverts = IntArray(self.nsnodes*self.numtractions)
        self.tractionvals = DoubleArray(constants.numberDegreesFreedom*self.numtractions)
        libpylith3d.read_tractions(
            self.tractionverts.ptr,    # intent(out)
            self.tractionvals.ptr,     # intent(out)
            &self.tscale,
            &self.numtractions,
            &self.nsnodes,
            &self.kr,
            self.tractionInputFile,
            &errorcode,
            errorstring,
            len(self.tractionInputFile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Read split node info
        self.nfault = IntArray(3*self.numfn)
        self.fault = DoubleArray(constants.numberDegreesFreedom*self.numfn)
        libpylith3d.read_split(
            self.fault.ptr,      # intent(out)
            self.nfault.ptr,     # intent(out)
            &self.numfn,
            &self.numflt,        # intent(out)
            &self.numnp,
            &self.numelv,
            &self.kr,
            self.spfile,
            &errorcode,
            errorstring,
            len(self.spfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Read slippery node info
        # Note that array Nslip is also required in functions sortmesh and sparsesetup
        # before it can be deallocated.
        self.nslip = IntArray(constants.numberSlipDimensions*self.numslp)
        libpylith3d.read_slip(
            self.nslip.ptr,      # intent(out)
            &self.numslp,
            &self.numsn,         # intent(out)
            &self.numnp,
            &self.iskopt,
            &self.kr,
            self.slfile,
            &errorcode,
            errorstring,
            len(self.slfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.idhist = IntArray(self.numnp)
        self.diforc = DoubleArray(constants.numberDegreesFreedom*self.numnp)
        libpylith3d.read_diff(
            self.diforc.ptr,    # intent(out)
            self.nslip.ptr,
            self.idhist.ptr,    # intent(out)
            &self.numslp,
            &self.numdif,
            &self.numnp,
            &self.kr,
            self.difile,
            &errorcode,
            errorstring,
            len(self.difile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)
        
        # Read Winkler forces and slippery Winkler forces.
        # All input is finished after this section.
        self.iwinkdef = IntArray(constants.numberDegreesFreedom*self.nwinke)
        self.iwinkid = IntArray(self.nwinke)
        self.winkdef = DoubleArray(constants.numberDegreesFreedom*self.nwinke)
        libpylith3d.read_wink(
            self.winkdef.ptr,     # intent(out)
            self._wscal,
            self.iwinkdef.ptr,    # intent(out)
            self.iwinkid.ptr,     # intent(out)
            &self.nwink,
            &self.nwinke,
            &self.kr,
            self.wfile,
            &errorcode,
            errorstring,
            len(self.wfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.iwinkxdef = IntArray(constants.numberDegreesFreedom*self.nwinkxe)
        self.iwinkxid = IntArray(self.nwinkxe)
        self.winkxdef = DoubleArray(constants.numberDegreesFreedom*self.nwinkxe)
        libpylith3d.read_wink(
            self.winkxdef.ptr,     # intent(out)
            self._wxscal,
            self.iwinkxdef.ptr,    # intent(out)
            self.iwinkxid.ptr,     # intent(out)
            &self.nwinkx,
            &self.nwinkxe,
            &self.kr,
            self.wxfile,
            &errorcode,
            errorstring,
            len(self.wxfile),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.trace.log("Hello from pl3dsetup.read (end)!")

        return


    cdef getdef(self):

        assert(self.ietypev <= elementtypes.numElementTypes)

        cdef int errorcode
        cdef char errorstring[errorstring_size]
        errorcode = 0
        
        cdef elementtypes.ElementType *elementType
        cdef int nsgauss

        elementType = &elementtypes.elementTypes[self.ietypev - 1]

        self.nen = elementType.nodes
        self.nsnodes = elementType.nodes2d

        if self.intord == 1:
            self.ngauss = elementType.fullGauss
            nsgauss = elementType.fullGauss2d
        elif self.intord == 2:
            self.ngauss = elementType.reducedGauss
            nsgauss = elementType.reducedGauss2d
        elif self.intord == 3:
            self.ngauss = elementType.fullGauss
            nsgauss = elementType.fullGauss2d
        else:
            raise ValueError("intord (%d) is not in [1, 2, 3]" % intord)

        cdef int nsd, nsd1
        nsd = constants.numberSpaceDimensions
        nsd1 = nsd + 1

        self.sh = DoubleArray(nsd1 * self.nen * self.ngauss)
        self.shj = DoubleArray(nsd1 * self.nen * self.ngauss)
        self.gauss = DoubleArray(nsd1 * self.ngauss)
            
        libpylith3d.preshape(
            self.sh.ptr,       # intent(out)
            self.shj.ptr,      # intent(out)
            self.gauss.ptr,    # intent(out)
            &self.intord,
            &self.ietypev,
            &self.nen,
            &self.ngauss,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.sh2d = DoubleArray(nsd * self.nsnodes * nsgauss)
        self.gauss2d = DoubleArray(nsd * nsgauss)

        libpylith3d.preshape2d(
            self.sh2d.ptr,       # intent(out)
            self.gauss2d.ptr,    # intent(out)
            &self.intord,
            &self.ietypev,
            &self.nsnodes,
            &nsgauss,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)
        
        self.connectivitySize = self.numelv*self.nen

        #
        # infetype array
        #
        
        self.nee = constants.numberDegreesFreedom * self.nen
        cdef int nec
        nec = 0

        cdef int *infetype
        infetype = &self.infetype[0] - 1 # convert to one-based
        infetype[1] = self.nen
        infetype[2] = self.ngauss
        infetype[3] = self.nee
        infetype[4] = nec
        assert(  4 == sizeof(self.infetype)/sizeof(self.infetype[0]))

        #
        # infetype2d array
        #

        cdef int numberSurfaceElementEquations, numberSurfaceElementCoordinates
        numberSurfaceElementEquations = (constants.numberDegreesFreedom *
                                         self.nsnodes)
        numberSurfaceElementCoordinates = (constants.numberSpaceDimensions *
                                           self.nsnodes)

        cdef int *infetype2d
        infetype2d = &self.infetype2d[0] - 1 # convert to one-based
        infetype2d[1] = self.nsnodes
        infetype2d[2] = nsgauss
        infetype2d[3] = numberSurfaceElementEquations # unused?
        infetype2d[4] = numberSurfaceElementCoordinates # unused?
        assert(    4 == sizeof(self.infetype2d)/sizeof(self.infetype2d[0]))

        return


    cdef numberequations(self):

        # This functions numbers equations based on BC and slippery node info.

        self.trace.log("Hello from pl3dsetup.numberequations (begin)!")
        
        print "Numbering global equations."

        # Create Idftn array for split nodes.  This can be deallocated
        # after meshwrite function has been called.
        self.idftn = IntArray(self.numflt)
        
        libpylith3d.id_split(
            self.nfault.ptr,
            self.idftn.ptr,     # intent(out)
            &self.numnp,
            &self.numfn,
            &self.numflt)

        # Determine global equations and store equation numbers in Id and Idx.
        self.id = IntArray(constants.numberSpaceDimensions*self.numnp)
        self.idx = IntArray(constants.numberSpaceDimensions*self.numnp)
        self.idslp = IntArray(self.numnp)

        # Number of equations
        
        libpylith3d.create_id(
            self.id.ptr,           # intent(out)
            self.idx.ptr,          # intent(out)
            self.ibond.ptr,
            self.nslip.ptr,
            self.idslp.ptr,        # intent(out)
            &self.numslp,
            &self.numnp,
            &self.numsn,
            &self.neq              # intent(out)
            )

        self.ipslp = IntArray(constants.numberSlipNeighbors*self.numsn)

        # If there are slippery nodes and the auto-rotation option is selected, find
        # neighboring nodes on the fault so that a best-fit plane can be determined at
        # each node.
        if self.numsn != 0 and self.iskopt == 2:
            self.nfind()

        # Assign appropriate equation numbers to Iwink array, and compact Wink
        # array to correspond to assigned BC.
        self.wink = DoubleArray(self.nwink)
        self.iwink = IntArray(2*self.nwink)

        libpylith3d.assign_wink(
            self.winkdef.ptr,
            self.wink.ptr,        # intent(out)
            self.iwinkdef.ptr,
            self.iwinkid.ptr,
            self.iwink.ptr,       # intent(out)
            self.id.ptr,
            &self.numnp,
            &self.nwink,
            &self.nwinke)

        # Assign appropriate equation numbers to Iwinkx array, and compact Winkx
        # array to correspond to assigned BC.
        self.winkx = DoubleArray(self.nwinkx)
        self.iwinkx = IntArray(2*self.nwinkx)

        libpylith3d.assign_wink(
            self.winkxdef.ptr,
            self.winkx.ptr,        # intent(out)
            self.iwinkxdef.ptr,
            self.iwinkxid.ptr,
            self.iwinkx.ptr,       # intent(out)
            self.idx.ptr,
            &self.numnp,
            &self.nwinkx,
            &self.nwinkxe)

        self.trace.log("Hello from pl3dsetup.numberequations (end)!")
            
        return


    cdef nfind(self):
        
        # Temporary arrays
        cdef DoubleArray xtmp
        cdef IntArray itmp, itmp1, itmp2
        
        xtmp = DoubleArray(self.numsn)
        itmp = IntArray(self.numsn)
        itmp1 = IntArray(self.numsn)
        itmp2 = IntArray(self.numsn)

        libpylith3d.nfind(
            self.x.ptr,
            xtmp.ptr,
            self.idslp.ptr,
            self.ipslp.ptr,    # intent(out)
            itmp.ptr,
            itmp1.ptr,
            itmp2.ptr,
            self.nslip.ptr,
            &self.numslp,
            &self.numsn,
            &self.numnp)

        return


    cdef sortmesh(self):

        # This function sorts elements into families and sorts all other items that are
        # affected by this.

        self.trace.log("Hello from pl3dsetup.sortmesh (begin)!")
        
        print "Renumbering elements, split nodes, and slippery nodes."

        self.sort_elements()

        # Sort split node entries.
        libpylith3d.sort_split_nodes(
            self.nfault.ptr,     # intent(inout)
            self.indxiel.ptr,
            &self.numfn,
            &self.numelv)

        # Sort slippery node entries.
        libpylith3d.sort_slip_nodes(
            self.nslip.ptr,      # intent(inout)
            self.indxiel.ptr,
            &self.numslp,
            &self.numelv)
            
        self.trace.log("Hello from pl3dsetup.sortmesh (end)!")

        return


    cdef sort_elements(self):
        # Sort elements into families.  The sorted elements are contained
        # in array Iens, and the index array for the new ordering is
        # Indxiel.  The index array for the original ordering is Ielindx.
        # The original element node array (Ien) and the associated
        # material type array (Mat) may be deallocated after sorting.
        
        self.iens = IntArray(self.nen*self.numelv)
        self.ivfamily = IntArray(6*self.nvfamilies)
        self.indxiel = IntArray(self.numelv)
        self.ielindx = IntArray(self.numelv)

        cdef IntArray ivftmp
        ivftmp = IntArray(self.nvfamilies)

        cdef int errorcode
        cdef char errorstring[errorstring_size]
        errorcode = 0
        
        libpylith3d.sort_elements(
            self.ien.ptr,
            self.mat.ptr,
            self.infmatmod,
            self.ivflist.ptr,
            self.ivfamily.ptr,    # intent(out)
            self.iens.ptr,        # intent(out)
            ivftmp.ptr,
            self.indxiel.ptr,     # intent(out)
            self.ielindx.ptr,     # intent(out)
            &self.nen,
            &self.ngauss,
            &self.maxvfamilies,
            &self.nvfamilies,
            &self.nprestrflag,
            &self.numelv,
            &self.numnp,
            &self.nstatesz,       # intent(out)
            &self.nstatesz0,      # intent(out)
            &self.npropsz,        # intent(out)
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)
        
        self.ien = None ### DEALLOC
        self.mat = None ### DEALLOC
        self.ivflist = None ### DEALLOC

        return


    cdef setupPETScLogging(self):
        
        self.autoprestrStage = PetscLogStageRegister("AutoPrestress Solve")
        self.elasticStage    = PetscLogStageRegister("Elastic Solve")
        self.viscousStage    = PetscLogStageRegister("Viscous Solve")
        
        self.iterateEvent = PetscLogEventRegister("Iterate", petsc.KSP_COOKIE)
        
        return


    cdef sparsesetup(self):

        # This function sets up sparse matrix and associated storage.

        cdef int errorcode
        cdef char errorstring[errorstring_size]
        errorcode = 0
        
        self.trace.log("Hello from pl3dsetup.sparsesetup (begin)!")
        
        print "Setting up sparse matrix storage."

        self.setupPETScLogging()

        # Arrays to map element equation numbers to global
        # Localize global equation numbers in element index arrays.

        self.lm = IntArray(constants.numberDegreesFreedom*self.connectivitySize)
        libpylith3d.local(
            self.id.ptr,
            &self.numnp,
            self.iens.ptr,
            self.lm.ptr,      # intent(out)
            &self.numelv,
            &self.nen)

        self.lmf = IntArray(self.connectivitySize)
        libpylith3d.localf(
            self.iens.ptr,
            self.lmf.ptr,      # intent(out)
            &self.numelv,
            self.nfault.ptr,
            &self.numfn,
            &self.nen)

        self.lmx = IntArray(constants.numberDegreesFreedom*self.connectivitySize)
        libpylith3d.localx(
            self.idx.ptr,
            &self.numnp,
            self.iens.ptr,
            self.lmx.ptr,      # intent(out)
            &self.numelv,
            self.nslip.ptr,
            &self.numslp,
            &self.nen)

        # Keeping this for now as it may be wanted for output
        if False: self.nslip = None ### DEALLOC

        # Allocate and populate sparse matrix arrays.  Some of these are
        # temporary and are then deleted after use.
        cdef int iwork
        iwork = 0
        libpylith3d.cmp_stiffsz(
            &self.neq,
            self.lm.ptr,
            self.lmx.ptr,
            &self.numelv,
            &iwork,          # intent(out)
            &self.numsn,
            &self.nen,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Temporary arrays
        cdef IntArray indx, link, nbrs
        indx = IntArray(self.neq)
        link = IntArray(iwork)
        nbrs = IntArray(iwork)

        cdef int nsizea

        libpylith3d.lnklst(
            &self.neq,
            self.lm.ptr,
            self.lmx.ptr,
            &self.numelv,
            &self.nen,
            &self.nee,
            indx.ptr,
            link.ptr,
            nbrs.ptr,
            &iwork,
            &nsizea,      # intent(out)
            &self.nnz,    # intent(out)
            &self.numsn,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self._mesh.createMat()

        libpylith3d.makemsr(
            &self._mesh.A,
            indx.ptr,
            link.ptr,
            nbrs.ptr,
            &self.neq,
            &self.nnz,
            &iwork,
            &self.nmin,
            &self.nmax,
            &self.wavg)

        print ""
        print ""
        print "Sparse matrix information:"
        print ""
        print "numberGlobalEquations:     %i" % self.neq
        print "workingArraySize:          %i" % iwork
        print "stiffnessMatrixSize:       %i" % (self.nnz-1)
        print "stiffnessOffDiagonalSize:  %i" % nsizea
        print "minimumNonzeroTermsPerRow: %i" % self.nmin
        print "maximumNonzeroTermsPerRow: %i" % self.nmax
        print "averageNonzeroTermsPerRow: %g" % self.wavg
        print ""
        
        self.trace.log("Hello from pl3dsetup.sparsesetup (end)!")

        return


    cdef allocateremaining(self):

        # This function allocates all remaining arrays that are needed for computations.
        
        self.trace.log("Hello from pl3dsetup.allocateremaining (begin)!")
        
        print "Allocating remaining storage."
        
        # Allocate memory for all additional arrays

        cdef int ntractflag, ngravflag, nconcflag, nextflag
        cdef int nwinkflag, nwinkxflag
        cdef int dim

        # Force vectors
        if self.numtractions != 0:
            ntractflag = 1
        else:
            ntractflag = 0
        
        ngravflag = 0
        dim = sizeof(self._grav)/sizeof(self._grav[0])
        for i from 0 <= i < dim:
            if self._grav[i] != 0.0:
                ngravflag = 1
                break

        if self.numberConcForces != 0 or self.numdif != 0:
            nconcflag = 1
        else:
            nconcflag = 0
        if ntractflag != 0 or ngravflag != 0 or nconcflag != 0:
            nextflag = 1
        else:
            nextflag = 0
        if self.nwink != 0:
            nwinkflag = 1
        else:
            nwinkflag = 0
        if self.nwinkx != 0:
            nwinkxflag = 1
        else:
            nwinkxflag = 0

        self.bextern = DoubleArray(nextflag*self.neq)
        self.btraction = DoubleArray(ntractflag*self.neq)
        self.bgravity = DoubleArray(ngravflag*self.neq)
        self.bconcForce = DoubleArray(nconcflag*self.neq)
        self.bwink = DoubleArray(nwinkflag*self.neq)
        self.bwinkx = DoubleArray(nwinkxflag*self.neq)
        self.bintern = DoubleArray(self.neq)
        self.bresid = DoubleArray(self.neq)
        self.dispVec = DoubleArray(self.neq)
        self.dprev = DoubleArray(self.nprevdflag*self.neq)
            
        # Displacement arrays
        self.d = DoubleArray(constants.numberDegreesFreedom*self.numnp)
        self.deld = DoubleArray(constants.numberDegreesFreedom*self.numnp)
        self.dcur = DoubleArray(constants.numberDegreesFreedom*self.numnp)

        # Slippery node arrays
        self.dx = DoubleArray(constants.numberDegreesFreedom*self.numnp)
        self.deldx = DoubleArray(constants.numberDegreesFreedom*self.numnp)
        self.dxcur = DoubleArray(constants.numberDegreesFreedom*self.numnp)

        # Split node arrays
        self.dfault = DoubleArray(constants.numberDegreesFreedom*self.numfn)
        self.tfault = DoubleArray(constants.numberDegreesFreedom*self.numfn)

        # Local stiffness matrix arrays
        self.s = DoubleArray(self.maxElementEquations*self.maxElementEquations)
        self.stemp = DoubleArray(self.maxElementEquations*self.maxElementEquations)

        # Element arrays
        self.state = DoubleArray(self.nstatesz)
        self.dstate = DoubleArray(self.nstatesz)
        self.dmat = DoubleArray(constants.materialMatrixDimension *
                                self.ngauss *
                                self.numelv)
        
        # This corresponds to the BLAS packed symmetric matrix format.
        listIddmat = [
             1, 2, 4, 7,11,16,
             2, 3, 5, 8,12,17,
             4, 5, 6, 9,13,18,
             7, 8, 9,10,14,19,
            11,12,13,14,15,20,
            16,17,18,19,20,21
            ]
        dim = sizeof(self.iddmat)/sizeof(self.iddmat[0])
        assert(len(listIddmat) == dim)
        for i from 0 <= i < dim:
            self.iddmat[i] = listIddmat[i]
        
        self.state0 = DoubleArray(self.nstatesz0)

        # Create arrays from lists that will be needed for the solution

        # nforce array ~ see nforce_def.inc
        cdef int *nforce
        nforce = &self.nforce[0] - 1 # convert to one-based
        nforce[1] = nextflag
        nforce[2] = ntractflag
        nforce[3] = ngravflag
        nforce[4] = nconcflag
        nforce[5] = self.nprestrflag
        nforce[6] = nwinkflag
        nforce[7] = nwinkxflag
        nforce[8] = self.nprevdflag
        assert(8 == sizeof(self.nforce)/sizeof(self.nforce[0]))
        
        # ncodat array ~ see ncodat_def.inc
        cdef int *ncodat
        ncodat = &self.ncodat[0] - 1 # convert to one-based
        ncodat[1] = self.icode
        ncodat[2] = self.idebug
        assert(2 == sizeof(self.ncodat)/sizeof(self.ncodat[0]))
            
        # npar array ~ see npar_def.inc
        cdef int *npar
        npar = &self.npar[0] - 1 # convert to one-based
        npar[1]  = self.numelv
        npar[2]  = self.numat
        npar[3]  = self.numtractions
        npar[4]  = self.numslp
        npar[5]  = self.numfn
        npar[6]  = self.ipstrs
        npar[7]  = self.ipauto
        npar[8]  = self.nstatesz
        npar[9]  = self.nstatesz0
        npar[10] = self.nvfamilies
        npar[11] = self.numdif
        npar[12] = self.intord
        assert(12 == sizeof(self.npar)/sizeof(self.npar[0]))

        # nprint array ~ see nprint_def.inc
        cdef int *nprint
        nprint = &self.nprint[0] - 1 # convert to one-based
        nprint[1] = self.icontr
        nprint[2] = self.idout
        nprint[3] = self.idsk
        nprint[4] = self.iucd
        assert(4 == sizeof(self.nprint)/sizeof(self.nprint[0]))

        # nsysdat array ~ see nsysdat_def.inc
        cdef int *nsysdat
        nsysdat = &self.nsysdat[0] - 1 # convert to one-based
        nsysdat[1]  = self.numnp
        nsysdat[2]  = self.neq
        nsysdat[3]  = self.nnz
        nsysdat[4]  = self.numrot
        nsysdat[5]  = self.nprestr
        nsysdat[6]  = self.numsn
        nsysdat[7]  = self.numflt
        nsysdat[8]  = self.npropsz
        nsysdat[9]  = self.nwink
        nsysdat[10] = self.nwinkx
        nsysdat[11] = self.iskopt
        assert( 11 == sizeof(self.nsysdat)/sizeof(self.nsysdat[0]))

        # nunits array ~ see nunits_def.inc
        cdef int *nunits
        nunits = &self.nunits[0] - 1 # convert to one-based
        nunits[1] = self.kti
        nunits[2] = self.kto
        nunits[3] = self.kr
        nunits[4] = self.kw
        nunits[5] = self.kp
        nunits[6] = self.kucd
        assert(6 == sizeof(self.nunits)/sizeof(self.nunits[0]))

        # nvisdat array ~ see nvisdat_def.inc
        cdef int *nvisdat
        nvisdat = &self.nvisdat[0] - 1 # convert to one-based
        nvisdat[1] = self.ncycle
        nvisdat[2] = self.nintg
        nvisdat[3] = self.lastep
        nvisdat[4] = self.nhist
        assert( 4 == sizeof(self.nvisdat)/sizeof(self.nvisdat[0]))
        
        # rgiter array ~ see rgiter_def.inc
        cdef double *rgiter
        rgiter = &self.rgiter[0] - 1 # convert to one-based
        rgiter[1] = self.stol
        rgiter[2] = self.dtol
        rgiter[3] = self.epert
        assert(3 == sizeof(self.rgiter)/sizeof(self.rgiter[0]))
        
        # rtimdat array ~ see rtimdat_def.inc
        cdef double *rtimdat
        rtimdat = &self.rtimdat[0] - 1 # convert to one-based
        cdef double deltp, alfap
        deltp = 0.0; alfap = 0.0
        rtimdat[1] = deltp
        rtimdat[2] = alfap
        rtimdat[3] = self.tpois
        rtimdat[4] = self.tyoungs
        assert( 4 == sizeof(self.rtimdat)/sizeof(self.rtimdat[0]))

        # ntimdat array ~ see ntimdat_def.inc
        cdef int *ntimdat
        cdef int nstep, maxitp, ntdinitp
        cdef int lgdefp, itmaxp, nittot
        cdef int nrftot, ndtot , ireform
        ntimdat = &self.ntimdat[0] - 1 # convert to one-based
        nstep  = 0; maxitp = 0; ntdinitp = 0
        lgdefp = 0; itmaxp = 0; nittot   = 0
        nrftot = 0; ndtot  = 0; ireform  = 0
        ntimdat[1] = nstep
        ntimdat[2] = maxitp
        ntimdat[3] = ntdinitp
        ntimdat[4] = lgdefp
        ntimdat[5] = itmaxp
        ntimdat[6] = nittot
        ntimdat[7] = nrftot
        ntimdat[8] = ndtot
        ntimdat[9] = ireform
        assert( 9 == sizeof(self.ntimdat)/sizeof(self.ntimdat[0]))

        self.trace.log("Hello from pl3dsetup.allocateremaining (end)!")

        return


    cdef meshwrite(self):

        # This function outputs mesh information.
        # In the near future, this needs to be broken into classes for
        # Ascii output, plot output, UCD output, etc.

        cdef int errorcode
        cdef char errorstring[errorstring_size]
        errorcode = 0

        self.trace.log("Hello from pl3dsetup.meshwriteascii (begin)!")
        
        print "Outputting Ascii mesh information."

        # Write out global parameters
        libpylith3d.write_global_info(
            self.title,
            &self.idout,
            &self.idsk,
            &self.numnp,
            &self.icode,
            &self.idebug,
            &self.kw,
            &self.kp,
            len(self.title))

        # Write out nodal coordinates
        libpylith3d.write_coords(
            self.x.ptr,
            &self.numnp,
            &self.kw,
            &self.kp,
            &self.idout,
            &self.idsk,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Write out nodal boundary condition info
        libpylith3d.write_bc(
            self.bond.ptr,
            self.ibond.ptr,
            &self.numnp,
            &self.kw,
            &self.idout,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Write out local coordinate rotations
        libpylith3d.write_skew(
            self.skew.ptr,
            &self.numrot,
            &self.iskopt,
            &self.numnp,
            &self.kw,
            &self.idout,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Write stress computation and subiteration parameters.
        libpylith3d.write_strscomp(
            &self.stol,
            &self.dtol,
            &self.epert,
            &self.kw,
            &self.idout)

        libpylith3d.write_subiter(
            &self.nprevdflag,
            &self.kw,
            &self.idout)

        # Write out time step information
        libpylith3d.write_timdat(
            self.delt.ptr,
            self.alfa.ptr,
            self.utol.ptr,
            self.ftol.ptr,
            self.etol.ptr,
            self.times.ptr,
            self.maxstp.ptr,
            self.maxit.ptr,
            self.ntdinit.ptr,
            self.lgdef.ptr,
            self.itmax.ptr,
            &self.nintg,
            &self.lastep,
            &self.kw,
            &self.idout,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Write out timesteps when full output is desired
        libpylith3d.write_fuldat(
            self.iprint.ptr,
            &self.icontr,
            &self.icode,
            &self.ncycle,
            &self.lastep,
            &self.kw,
            &self.kp,
            &self.idout,
            &self.idsk,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Write out state variables desired for output
        libpylith3d.write_stateout(
            self.istatout,
            self.nstatout,
            &self.kw,
            &self.kp,
            &self.idout,
            &self.idsk,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Write out load history information and deallocate 'times' array
        libpylith3d.write_hist(
            self.histry.ptr,
            self.times.ptr,
            &self.nhist,
            &self.lastep,
            &self.kw,
            &self.idout,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.times = None ### DEALLOC

        # Write element info
        libpylith3d.write_element_info(
            &self.numelv,
            &self.nen,
            &self.ngauss,
            &self.ietypev,
            &self.intord,
            &self.ipstrs,
            &self.ipauto,
            &self.tpois,
            &self.tyoungs,
            &self.kw,
            &self.idout)

        # Write element node array and deallocate 'indxiel'
        libpylith3d.write_connect(
            self.iens.ptr,
            self.ivfamily.ptr,
            self.indxiel.ptr,
            &self.nen,
            &self.ngauss,
            &self.numelv,
            &self.ietypev,
            &self.nvfamilies,
            &self.kw,
            &self.kp,
            &self.idout,
            &self.idsk,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.indxiel = None ### DEALLOC

        # Write material properties
        libpylith3d.write_props(
            self._prop.ptr,
            self._grav,
            self.ivfamily.ptr,
            self.infmatmod,
            &self.nvfamilies,
            &self.npropsz,
            &self.idout,
            &self.idsk,
            &self.kw,
            &self.kp,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Write mesh info to UCD file, if requested
        if self.iucd >= 0:
            libpylith3d.write_ucd_mesh(
                self.x.ptr,
                &self.numnp,
                self.iens.ptr,
                self.ivfamily.ptr,
                &self.numelv,
                &self.nvfamilies,
                self.sh.ptr,
                &self.nen,
                &self.ngauss,
                &self.ietypev,
                self.istatout,
                self.nstatout,
                &self.kucd,
                &self.iucd,
                self.ucdroot,
                len(self.ucdroot))

        # Write traction info
        libpylith3d.write_tractions(
            self.tractionverts.ptr,
            self.tractionvals.ptr,
            &self.numtractions,
            &self.nsnodes,
            &self.kw,
            &self.idout,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)
   
        # Write split node info
        libpylith3d.write_split(
            self.fault.ptr,
            self.nfault.ptr,
            &self.numfn,
            &self.kw,
            &self.kp,
            &self.idout,
            &self.idsk,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Write slippery node info
        libpylith3d.write_slip(
            self.nslip.ptr,
            &self.numslp,
            &self.numsn,
            &self.kw,
            &self.kp,
            &self.idout,
            &self.idsk,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        # Write differential force info and deallocate 'nslip'
        libpylith3d.write_diff(
            self.diforc.ptr,
            self.nslip.ptr,
            self.idhist.ptr,
            &self.numslp,
            &self.numdif,
            &self.numnp,
            &self.kw,
            &self.idout,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.nslip = None ### DEALLOC

        # Write split nodes to plot file, if requested and deallocate 'idftn'
        libpylith3d.write_split_plot(
            self.idftn.ptr,
            &self.numflt,
            &self.kp,
            &self.idsk)

        self.idftn = None ### DEALLOC

        # Write Winkler force info and deallocate definition arrays
        libpylith3d.write_wink(
            self.winkdef.ptr,
            self.iwinkdef.ptr,
            self.iwinkid.ptr,
            &self.nwinke,
            &self.kw,
            &self.idout,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.winkdef = None ### DEALLOC
        self.iwinkdef = None ### DEALLOC

        # Write slippery node Winkler force info and deallocate definition arrays
        libpylith3d.write_winkx(
            self.winkxdef.ptr,
            self.iwinkxdef.ptr,
            self.iwinkxid.ptr,
            &self.nwinkxe,
            &self.kw,
            &self.idout,
            &errorcode,
            errorstring,
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)

        self.winkxdef = None ### DEALLOC
        self.iwinkxdef = None ### DEALLOC

        # Write sparse matrix info
        libpylith3d.write_sparse_info(
            &self.neq,
            &self.nnz,
            &self.nmin,
            &self.nmax,
            &self.wavg,
            &self.idout,
            &self.kw)

        self.trace.log("Hello from pl3dsetup.meshwrite (end)!")

        return


# The function of this code is to call the elastic and time-dependent solution
# drivers.  To do this, a number of previously-defined parameters need to be
# bundled into lists.


    cdef solveElastic(self):
        
        cdef int errorcode
        cdef char errorstring[errorstring_size]
        errorcode = 0
        
        libpylith3d.elastc(
            &self._mesh.A,         # sparse
            &self._mesh.rhs,
            &self._mesh.sol,
            self.bextern.ptr,      # force
            self.btraction.ptr,
            self.bgravity.ptr,
            self.bconcForce.ptr,
            self.bintern.ptr,
            self.bresid.ptr,
            self.bwink.ptr,
            self.bwinkx.ptr,
            self.dispVec.ptr,
            self.dprev.ptr,
            self.nforce,
            self._grav,
            self.x.ptr,            # global
            self.d.ptr,
            self.deld.ptr,
            self.dcur.ptr,
            self.id.ptr,
            self.iwink.ptr,
            self.wink.ptr,
            self.nsysdat,
            self.iddmat,
            self.ibond.ptr,        # BC
            self.bond.ptr,
            self.dx.ptr,           # slip
            self.deldx.ptr,
            self.dxcur.ptr,
            self.diforc.ptr,
            self.idx.ptr,
            self.iwinkx.ptr,
            self.winkx.ptr,
            self.idslp.ptr,
            self.ipslp.ptr,
            self.idhist.ptr,
            self.fault.ptr,        # fault
            self.nfault.ptr,
            self.dfault.ptr,
            self.tfault.ptr,
            self.s.ptr,            # stiff
            self.stemp.ptr,
            self.state.ptr,        # element
            self.dstate.ptr,
            self.state0.ptr,
            self.dmat.ptr,
            self.iens.ptr,
            self.lm.ptr,
            self.lmx.ptr,
            self.lmf.ptr,
            self.ivfamily.ptr,
            self.npar,
            self.ielindx.ptr,
            self.tractionverts.ptr, # traction
            self.tractionvals.ptr,
            self.gauss2d.ptr,
            self.sh2d.ptr,
            self.infetype2d,
            self._prop.ptr,        # material
            self.infmatmod,
            self.gauss.ptr,        # eltype
            self.sh.ptr,
            self.shj.ptr,
            self.infetype,
            self.histry.ptr,       # timdat
            self.rtimdat,
            self.ntimdat,
            self.nvisdat,
            self.maxstp.ptr,
            self.delt.ptr,
            self.alfa.ptr,
            self.maxit.ptr,
            self.ntdinit.ptr,
            self.lgdef.ptr,
            self.utol.ptr,
            self.ftol.ptr,
            self.etol.ptr,
            self.itmax.ptr,
            self.rgiter,           # stresscmp
            self.skew.ptr,         # skew
            self.ncodat,           # ioinfo
            self.nunits,
            self.nprint,
            self.istatout,
            self.nstatout,
            self.ucdroot,          # files
            &self.elasticStage,    # PETSc logging
            &self.iterateEvent,
            &errorcode,
            errorstring,
            len(self.ucdroot),
            sizeof(errorstring))
        exceptionhandler(errorcode, errorstring)
            
        return

    cdef interpolatePoints(self, points):
        return self._mesh.interpolatePoints(points)

    cdef viscos_setup(self):
        return

    cdef runSimulation(self):
        # First define all of the lists that maintain variable values.  The
        # variables in these lists are altered during the running of the code
        # and should not be accessed directly except as a member of the list.
        # They should not have been defined previously.

        cdef int errorcode
        cdef char errorstring[errorstring_size]
        errorcode = 0
        
        self.trace.log("Hello from pl3drun.run (begin)!")
        
        print "Beginning problem solution."

        if False: # Temporarily out-of-order
            # Output approximate memory usage
            self.memorySizeMB =0.0
            self.memorySizeMB=self.memorySize/(1024.0*1024.0)

            print ""
            print "Approximate memory allocation for f77 arrays (MB): %g" % self.memorySizeMB

        # Compute gravitational prestresses, if requeste
        if self.icode == 2 or self.icode == 3:  # elasticSolution or fullSolution
            if self.ipstrs == 1:
                libpylith3d.autoprestr(
                    &self._mesh.A,         # sparse
                    &self._mesh.rhs,
                    &self._mesh.sol,
                    self.bextern.ptr,      # force
                    self.btraction.ptr,
                    self.bgravity.ptr,
                    self.bconcForce.ptr,
                    self.bintern.ptr,
                    self.bresid.ptr,
                    self.bwink.ptr,
                    self.bwinkx.ptr,
                    self.dispVec.ptr,
                    self.dprev.ptr,
                    self.nforce,
                    self._grav,
                    self.x.ptr,            # global
                    self.d.ptr,
                    self.deld.ptr,
                    self.dcur.ptr,
                    self.id.ptr,
                    self.iwink.ptr,
                    self.wink.ptr,
                    self.nsysdat,
                    self.iddmat,
                    self.ibond.ptr,        # BC
                    self.bond.ptr,
                    self.dx.ptr,           # slip
                    self.deldx.ptr,
                    self.dxcur.ptr,
                    self.diforc.ptr,
                    self.idx.ptr,
                    self.iwinkx.ptr,
                    self.winkx.ptr,
                    self.idslp.ptr,
                    self.ipslp.ptr,
                    self.idhist.ptr,
                    self.fault.ptr,        # split
                    self.nfault.ptr,
                    self.dfault.ptr,
                    self.tfault.ptr,
                    self.s.ptr,            # stiff
                    self.stemp.ptr,
                    self.state.ptr,        # element
                    self.dstate.ptr,
                    self.state0.ptr,
                    self.dmat.ptr,
                    self.iens.ptr,
                    self.lm.ptr,
                    self.lmx.ptr,
                    self.lmf.ptr,
                    self.ivfamily.ptr,
                    self.npar,
                    self.ielindx.ptr,
                    self.tractionverts.ptr, # traction
                    self.tractionvals.ptr,
                    self.gauss2d.ptr,
                    self.sh2d.ptr,
                    self.infetype2d,
                    self._prop.ptr,        # material
                    self.infmatmod,
                    self.gauss.ptr,        # eltype
                    self.sh.ptr,
                    self.shj.ptr,
                    self.infetype,
                    self.histry.ptr,       # timdat
                    self.rtimdat,
                    self.ntimdat,
                    self.nvisdat,
                    self.maxstp.ptr,
                    self.delt.ptr,
                    self.alfa.ptr,
                    self.maxit.ptr,
                    self.ntdinit.ptr,
                    self.lgdef.ptr,
                    self.utol.ptr,
                    self.ftol.ptr,
                    self.etol.ptr,
                    self.itmax.ptr,
                    self.rgiter,           # stresscmp
                    self.skew.ptr,         # skew
                    self.ncodat,           # ioinfo
                    self.nunits,
                    self.nprint,
                    self.istatout,
                    self.nstatout,
                    self.ucdroot,          # files
                    &self.autoprestrStage, # PETSc logging
                    &self.iterateEvent,
                    &errorcode,
                    errorstring,
                    len(self.ucdroot),
                    sizeof(errorstring))
                exceptionhandler(errorcode, errorstring)

            # Perform elastic solution, if requested.
            self.solveElastic()
            self._mesh.outputMesh(self.fileRoot)

        # Perform time-dependent solution, if requested.

        cdef int numCycles, numTimeStepGroups, numslp
        cdef int iskopt, icontr, indexx, totalSteps
        
        cdef int nextStartStep, timeStep, startStep
        cdef double time

        cdef int maxitp, ntdinitp, lgdefp, itmaxp
        cdef double dt, alfap, gtol[3]

        cdef int ltim, cycle, tsGroup, j, skc
        
        if self.icode == 3 and self.nintg > 1:
            if self.pythonTimestep:
                PetscLogStagePush(self.viscousStage)
                
                numCycles         = self.nvisdat[0]
                numTimeStepGroups = self.nvisdat[1]
                numslp            = self.npar[3]
                iskopt            = self.nsysdat[10]
                icontr            = self.nprint[0]
                indexx            = 1 # Fortran index
                totalSteps        = 0 # This is ntot
                
                for cycle from 0 <= cycle < numCycles:
                    
                    if numCycles > 1: print '     working on cycle %d' % cycle
                    nextStartStep = 0 # This is naxstp
                    timeStep      = 0 # This is nstep
                    startStep     = 0 # This is nfirst
                    time          = 0.0

                    for tsGroup from 1 <= tsGroup < numTimeStepGroups:
                        # Define constants
                        dt = self.delt.ptr[tsGroup] # This is deltp
                        self.rtimdat[0] = dt
                        alfap = self.alfa.ptr[tsGroup]
                        self.rtimdat[1] = alfap
                        self.ntimdat[0] = timeStep
                        maxitp = self.maxit.ptr[tsGroup]
                        self.ntimdat[1] = maxitp
                        ntdinitp = self.ntdinit.ptr[tsGroup]
                        self.ntimdat[2] = ntdinitp
                        lgdefp = self.lgdef.ptr[tsGroup]
                        self.ntimdat[3] = lgdefp
                        itmaxp = self.itmax.ptr[tsGroup]
                        self.ntimdat[4] = itmaxp
                        gtol[0] = self.utol.ptr[tsGroup]
                        gtol[1] = self.ftol.ptr[tsGroup]
                        gtol[2] = self.etol.ptr[tsGroup]
                        startStep     = nextStartStep + 1
                        nextStartStep = startStep + self.maxstp.ptr[tsGroup] - 1

                        ltim = 1

                        for j from startStep <= j < nextStartStep+1:
                            totalSteps = totalSteps + 1
                            timeStep   = timeStep  + 1
                            self.ntimdat[0] = timeStep
                            time = time + dt
                            skc   = (numslp != 0 and (iskopt == 2 or (iskopt <= 0 and abs(iskopt) == timeStep)))

                            libpylith3d.viscos_step(
                                &self._mesh.A,         # sparse
                                &self._mesh.rhs,
                                &self._mesh.sol,
                                self.bextern.ptr,      # force
                                self.btraction.ptr,
                                self.bgravity.ptr,
                                self.bconcForce.ptr,
                                self.bintern.ptr,
                                self.bresid.ptr,
                                self.bwink.ptr,
                                self.bwinkx.ptr,
                                self.dispVec.ptr,
                                self.dprev.ptr,
                                self.nforce,
                                self._grav,
                                self.x.ptr,            # global
                                self.d.ptr,
                                self.deld.ptr,
                                self.dcur.ptr,
                                self.id.ptr,
                                self.iwink.ptr,
                                self.wink.ptr,
                                self.nsysdat,
                                self.iddmat,
                                self.ibond.ptr,        # BC
                                self.bond.ptr,
                                self.dx.ptr,           # slip
                                self.deldx.ptr,
                                self.dxcur.ptr,
                                self.diforc.ptr,
                                self.idx.ptr,
                                self.iwinkx.ptr,
                                self.winkx.ptr,
                                self.idslp.ptr,
                                self.ipslp.ptr,
                                self.idhist.ptr,
                                self.fault.ptr,        # fault
                                self.nfault.ptr,
                                self.dfault.ptr,
                                self.tfault.ptr,
                                self.s.ptr,            # stiff
                                self.stemp.ptr,
                                self.state.ptr,        # element
                                self.dstate.ptr,
                                self.state0.ptr,
                                self.dmat.ptr,
                                self.iens.ptr,
                                self.lm.ptr,
                                self.lmx.ptr,
                                self.lmf.ptr,
                                self.ivfamily.ptr,
                                self.npar,
                                self.ielindx.ptr,
                                self.tractionverts.ptr, # traction
                                self.tractionvals.ptr,
                                self.gauss2d.ptr,
                                self.sh2d.ptr,
                                self.infetype2d,
                                self._prop.ptr,        # material
                                self.infmatmod,
                                self.gauss.ptr,        # eltype
                                self.sh.ptr,
                                self.shj.ptr,
                                self.infetype,
                                self.histry.ptr,       # timdat
                                self.rtimdat,
                                self.ntimdat,
                                self.nvisdat,
                                self.maxstp.ptr,
                                self.delt.ptr,
                                self.alfa.ptr,
                                self.maxit.ptr,
                                self.ntdinit.ptr,
                                self.lgdef.ptr,
                                self.utol.ptr,
                                self.ftol.ptr,
                                self.etol.ptr,
                                self.itmax.ptr,
                                self.rgiter,           # stresscmp
                                self.skew.ptr,         # skew
                                self.iprint.ptr,       # ioinfo
                                self.ncodat,
                                self.nunits,
                                self.nprint,
                                self.istatout,
                                self.nstatout,
                                self.ucdroot,          # files
                                &self.viscousStage,    # PETSc logging
                                &self.iterateEvent,
                                &totalSteps,
                                &ltim,
                                &indexx,
                                &cycle,
                                &tsGroup,
                                &j,
                                &skc,
                                &startStep,
                                &timeStep,
                                &time,
                                &dt,
                                &lgdefp,
                                gtol,
                                &errorcode,
                                errorstring,
                                len(self.ucdroot),
                                sizeof(errorstring))
                            exceptionhandler(errorcode, errorstring)
                            
                            ltim = 0
                            if (totalSteps == self.iprint.ptr[indexx-1]):
                                self._mesh.outputMesh(self.fileRoot+'.'+str(totalSteps))
                                indexx = indexx +  1
                            if (indexx > icontr): indexx = icontr

                print " Total number of equilibrium iterations        =",self.ntimdat[5]
                print " Total number of stiffness matrix reformations =",self.ntimdat[6]
                print " Total number of displacement subiterations    =",self.ntimdat[7]
                
                libpylith3d.viscos_cleanup(
                    self.ntimdat,
                    self.nprint,
                    self.nunits,
                    &errorcode,
                    errorstring,
                    sizeof(errorstring))
                exceptionhandler(errorcode, errorstring)
                
                PetscLogStagePop()
                
            else:
                
                libpylith3d.viscos(
                    &self._mesh.A,         # sparse
                    &self._mesh.rhs,
                    &self._mesh.sol,
                    self.bextern.ptr,      # force
                    self.btraction.ptr,
                    self.bgravity.ptr,
                    self.bconcForce.ptr,
                    self.bintern.ptr,
                    self.bresid.ptr,
                    self.bwink.ptr,
                    self.bwinkx.ptr,
                    self.dispVec.ptr,
                    self.dprev.ptr,
                    self.nforce,
                    self._grav,
                    self.x.ptr,            # global
                    self.d.ptr,
                    self.deld.ptr,
                    self.dcur.ptr,
                    self.id.ptr,
                    self.iwink.ptr,
                    self.wink.ptr,
                    self.nsysdat,
                    self.iddmat,
                    self.ibond.ptr,        # BC
                    self.bond.ptr,
                    self.dx.ptr,           # slip
                    self.deldx.ptr,
                    self.dxcur.ptr,
                    self.diforc.ptr,
                    self.idx.ptr,
                    self.iwinkx.ptr,
                    self.winkx.ptr,
                    self.idslp.ptr,
                    self.ipslp.ptr,
                    self.idhist.ptr,
                    self.fault.ptr,        # fault
                    self.nfault.ptr,
                    self.dfault.ptr,
                    self.tfault.ptr,
                    self.s.ptr,            # stiff
                    self.stemp.ptr,
                    self.state.ptr,        # element
                    self.dstate.ptr,
                    self.state0.ptr,
                    self.dmat.ptr,
                    self.iens.ptr,
                    self.lm.ptr,
                    self.lmx.ptr,
                    self.lmf.ptr,
                    self.ivfamily.ptr,
                    self.npar,
                    self.ielindx.ptr,
                    self.tractionverts.ptr, # traction
                    self.tractionvals.ptr,
                    self.gauss2d.ptr,
                    self.sh2d.ptr,
                    self.infetype2d,
                    self._prop.ptr,        # material
                    self.infmatmod,
                    self.gauss.ptr,        # eltype
                    self.sh.ptr,
                    self.shj.ptr,
                    self.infetype,
                    self.histry.ptr,       # timdat
                    self.rtimdat,
                    self.ntimdat,
                    self.nvisdat,
                    self.maxstp.ptr,
                    self.delt.ptr,
                    self.alfa.ptr,
                    self.maxit.ptr,
                    self.ntdinit.ptr,
                    self.lgdef.ptr,
                    self.utol.ptr,
                    self.ftol.ptr,
                    self.etol.ptr,
                    self.itmax.ptr,
                    self.rgiter,           # stresscmp
                    self.skew.ptr,         # skew
                    self.iprint.ptr,       # ioinfo
                    self.ncodat,
                    self.nunits,
                    self.nprint,
                    self.istatout,
                    self.nstatout,
                    self.ucdroot,          # files
                    &self.viscousStage,    # PETSc logging
                    &self.iterateEvent,
                    &errorcode,
                    errorstring,
                    len(self.ucdroot),
                    sizeof(errorstring))
                exceptionhandler(errorcode, errorstring)

        self._mesh.destroyMat()

        self.trace.log("Hello from pl3drun.run (end)!")
        
        return



def try_binio(kucd):
    cdef int u
    cdef int errorcode
    cdef char errorstring[errorstring_size]

    u = kucd
    errorcode = 0
    
    libpylith3d.try_binio(&u, &errorcode, errorstring, sizeof(errorstring))
    exceptionhandler(errorcode, errorstring)
    
    return



# end of file 
