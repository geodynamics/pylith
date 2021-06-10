# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/apps/PyLithApp.py
##
# @brief Python PyLith application

from .PetscApplication import PetscApplication

# PyLithApp class

def problemFactory(name):
    """Factory for problem items.
    """
    from pythia.pyre.inventory import facility
    from pylith.problems.TimeDependent import TimeDependent
    return facility(name, family="problem", factory=TimeDependent)



class PyLithApp(PetscApplication):
    """Python PyLithApp application.
    """

    import pythia.pyre.inventory

    pdbOn = pythia.pyre.inventory.bool("start_python_debugger", default=False)
    pdbOn.meta['tip'] = "Start python debugger at beginning of main()."

    typos = pythia.pyre.inventory.str("typos", default="pedantic",
                                      validator=pythia.pyre.inventory.choice(['relaxed', 'strict', 'pedantic']))
    typos.meta['tip'] = "Specifies the handling of unknown properties and facilities"

    from pylith.utils.SimulationMetadata import SimulationMetadata
    metadata = pythia.pyre.inventory.facility(
        "metadata", family="simulation_metadata", factory=SimulationMetadata)
    metadata.meta["tip"] = "Simulation metadata."

    initializeOnly = pythia.pyre.inventory.bool(
        "initialize_only", default=False)
    initializeOnly.meta['tip'] = "Stop simulation after initializing problem."

    from pylith.utils.DumpParametersJson import DumpParametersJson
    parameters = pythia.pyre.inventory.facility(
        "dump_parameters", family="dump_parameters", factory=DumpParametersJson)
    parameters.meta['tip'] = "Dump parameters used and version information to file."

    from pylith.problems.SingleProblem import SingleProblem
    problems = pythia.pyre.inventory.facilityArray(
        "problems", itemFactory=problemFactory, factory=SingleProblem)
    problems.meta['tip'] = "Computational problem(s) to solve."
    
    # initialize 1 problem (original text)
    # problem = pythia.pyre.inventory.facility(
    #     "problem", family="problem", factory=TimeDependent)
    #problem.meta['tip'] = "Computational problem to solve."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="pylithapp"):
        """Constructor.
        """
        PetscApplication.__init__(self, name)
        self._loggingPrefix = "PyLith "
        return

    def main(self, *args, **kwds):
        """Run the application.
        """
        if self.pdbOn:
            import pdb
            pdb.set_trace()

        # Dump parameters and version information
        self.parameters.preinitialize()
        self.parameters.write(self)

        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Running on %d process(es)." % comm.size)

        from pylith.utils.profiling import resourceUsageString
        self._debug.log(resourceUsageString())

        self._setupLogging()


        # Create mesh for each problem
        self._eventLogger.stagePush("Meshing")
        for problem in self.problems.components():
            problem.createMesh()
        self._debug.log(resourceUsageString())
        self._eventLogger.stagePop()

        # Setup problem, verify configuration, and then initialize
        self._eventLogger.stagePush("Setup")
        for problem in self.problems.components():
            problem.preinitialize()
            self._debug.log(resourceUsageString())

            problem.verifyConfiguration()

            problem.initialize()
            self._debug.log(resourceUsageString())
        self._eventLogger.stagePop()

        # If initializing only, stop before running problem
        if self.initializeOnly:
            return

        # Cycle over all problems
        # TODO: load numIterations from input file
        # TODO: add control over integration time (start and end times)
        numIterations = 1 # number of earthquake cycles !!! later import this from an input file
        count = 0
        while count < numIterations:
            for problem in self.problems.components():
                problem.run()
                self._debug.log(resourceUsageString())
            count += 1

        # Cleanup
        self._eventLogger.stagePush("Finalize")
        for problem in self.problems.components():
            problem.finalize()
        self._eventLogger.stagePop()

        return

    def version(self):
        from pylith.utils.CollectVersionInfo import CollectVersionInfo
        msg = CollectVersionInfo.asString()
        msg += "\n"

        # Citation information
        msg += "If you publish results based on computations with PyLith please cite the following:\n" \
            "(use --include-citations during your simulations to display a list specific to your computation):\n\n"
        for citation in self.citations():
            msg += citation + "\n"

        print(msg)
        return

    def citations(self):
        import pylith.utils.utils as utils
        v = utils.PylithVersion()
        verNum = v.version()
        verYear = 2017
        verDOI = v.doi()

        software = ("@Manual{PyLith:software,\n"
                    "  title        = {PyLith v%s},\n"
                    "  author       = {Aagaard, B. and Knepley, M. and Williams, C.},\n"
                    "  organization = {Computational Infrastructure for Geodynamics (CIG)},\n"
                    "  address      = {University of California, Davis},\n"
                    "  year         = {%d},\n"
                    "  doi         = {http://doi.org/%s}\n"
                    "}\n" % (verNum, verYear, verDOI)
                    )

        manual = ("@Manual{PyLith:manual,\n"
                  "  title        = {PyLith User Manual, Version %s},\n"
                  "  author       = {Aagaard, B. and Knepley, M. and Williams, C.},\n"
                  "  organization = {Computational Infrastructure for Geodynamics (CIG)},\n"
                  "  address      = {University of California, Davis},\n"
                  "  year         = {%d},\n"
                  "  note         = {http://www.geodynamics.org/cig/software/pylith/pylith\_manual-%s.pdf}\n"
                  "}\n" % (verNum, verYear, verNum)
                  )

        faultRup = ("@Article{Aagaard:Knepley:Williams:JGR:2013,\n"
                    "  author   = {Aagaard, B.~T. and Knepley, M.~G. and Wiliams, C.~A.},\n"
                    "  title    = {A domain decomposition approach to implementing fault slip "
                    "in finite-element models of quasi-static and dynamic crustal deformation},\n"
                    "  journal  = {Journal of Geophysical Research Solid Earth},\n"
                    "  year     = {2013},\n"
                    "  volume   = {118},\n"
                    "  pages    = {3059--3079},\n"
                    "  doi      = {10.1002/jgrb.50217}\n"
                    "}\n"
                    )

        entries = (software, manual, faultRup)
        return entries

    def showHelp(self):
        msg = (
            "Before you ask for help, consult the PyLith user manual and try to debug on your own.\n"
            "You will likely find other useful information while making progress on your original issue.\n"
            "\n"
            "Helpful Resources:\n"
            "* User manual (https://geodynamics.org/cig/software/pylith/)\n"
            "* PyLith Tutorials (https://wiki.geodynamics.org/software:pylith:start)\n"
            "* pylithinfo script\n"
            "    Running pylithinfo --verbose [-o pylith_parameters.txt] [PyLith args]\n"
            "    will dump all parameters with descriptions to pylith_parameters.txt.\n"
            "\n"
            "Add $PYLITH_DIR/share/settings/petsc_monitor.cfg to your command line arguments\n"
            "to turn on several PETSc monitors:"
            "  pylith YOUR_FILE.cfg PATH_TO_PYITH_SHARE/share/settings/petsc_monitor.cfg\n"
            "\n"
            "If you still need help, visit the PyLith category on the CIG community forum:\n"
            "https://community.geodynamics.org.\n"
            "When submitting a question about running a simulation, be sure to include the info:\n"
            "\n"
            "1. Describe what you are trying to do\n"
            "  a. Overview of the problem and boundary conditions (diagrams are very helpful)\n"
            "  b. 2-D or 3-D\n"
            "  c. Cell type (tri, quad, hex, or tet)\n"
            "  d. Type of fault: prescribed slip or spontaneous rupture\n"
            "2. Attach the PyLith parameters .json file\n"
            "3. Send the *entire* error message, not just what you think is important (entire log is best).\n"
            "\n"
            "Description and help for PyLithApp component:\n"
        )
        if self.inventory.usage:
            print(msg)

        PetscApplication.showHelp(self)

        msg = (
            "\nExamples using step01.cfg in directory examples/3d/hex8):\n"
            "1. List components and properties for a given component (--help)\n"
            "  pylith step01.cfg --problem.bc.z_neg.help\n"
            "\n"
            "2. List components of a given component (--help-components)\n"
            "  pylith step01.cfg --problem.bc.z_neg.help-components\n"
            "\n"
            "3. List properties of a given component (--help-properties)\n"
            "  pylith step01.cfg --problem.bc.z_neg.help-properties\n"
        )
        if self.inventory.usage:
            print(msg)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        PetscApplication._configure(self)
        return

    def _setupLogging(self):
        """Setup event logging.
        """
        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("PyLith")
        logger.initialize()

        self._eventLogger = logger
        return


# ======================================================================
# Local version of InfoApp that only configures itself. Workaround for
# not adding --help-all (or similar property) in Pyre Application.
#
# Note: We want --help-all to display settings before launching a
# parallel job.
class InfoApp(PyLithApp):

    def __init__(self, args, name="pylithapp"):
        """Constructor.
        """
        PyLithApp.__init__(self, name)
        self.pylithargs = args
        return

    def onLoginNode(self, *args, **kwds):
        """Instead of scheduling job, do nothing.
        """
        return

    def getArgv(self, *args, **kwds):
        """Prevent PyLith from getting all of the command line arguments. Use
        only the ones relevant to PyLith which are specified in the arg to
        the constructor.
        """
        argv = kwds.get('argv')
        if argv is None:
            argv = self.pylithargs
        else:
            self.arg0 = argv[0]
            self._requires = kwds.get('requires')
            argv = argv[1:]
        return argv


# End of file
