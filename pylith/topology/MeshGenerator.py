# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from pylith.utils.PetscComponent import PetscComponent


class MeshGenerator(PetscComponent):
    """
    Abstract base class for mesh generator.
    """

    import pythia.pyre.inventory

    debug = pythia.pyre.inventory.bool("debug", default=False)
    debug.meta['tip'] = "Debugging flag for mesh."

    interpolate = pythia.pyre.inventory.bool("interpolate", default=True)
    interpolate.meta['tip'] = "Build intermediate mesh topology elements"

    def __init__(self, name="meshgenerator"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="meshgenerator")
        self.debug = False
        self.interpolate = True
        return

    def preinitialize(self, problem):
        """Do minimal initialization.
        """

    def create(self, normalizer, faults=None):
        """Generate a Mesh.
        """

        # Need to nondimensionalize coordinates.

        raise NotImplementedError("MeshGenerator.create() not implemented.")

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)

    def _adjustTopology(self, mesh, interfaces, problem):
        """Adjust topology for interface implementation.
        """
        logEvent = "%sadjTopo" % self._loggingPrefix
        self._eventLogger.eventBegin(logEvent)

        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()

        if not interfaces is None:
            for interface in interfaces:
                if 0 == comm.rank:
                    self._info.log("Adjusting topology for fault '%s'." % interface.labelName)
                interface.preinitialize(problem)
                interface.adjustTopology(mesh)

        self._eventLogger.eventEnd(logEvent)

    def _setupLogging(self):
        """Setup event logging.
        """
        if not "_loggingPrefix" in dir(self):
            self._loggingPrefix = ""

        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("Mesh Generator")
        logger.initialize()

        events = ["create",
                  "adjTopo"]
        for event in events:
            logger.registerEvent("%s%s" % (self._loggingPrefix, event))

        self._eventLogger = logger


# End of file
