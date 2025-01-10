# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.utils.PetscComponent import PetscComponent


class MeshGenerator(PetscComponent):
    """
    Abstract base class for mesh generator.
    """

    def __init__(self, name="meshgenerator"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="meshgenerator")

    def preinitialize(self, problem):
        """Do minimal initialization.
        """

    def create(self, normalizer, faults=None):
        """Generate a Mesh.
        """

        # Need to nondimensionalize coordinates.

        raise NotImplementedError("MeshGenerator.create() not implemented.")    

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)

    def _adjustTopology(self, mesh, interfaces, problem):
        """Adjust topology for interface implementation.
        """
        logEvent = f"{self._loggingPrefix}adjTopo"
        self._eventLogger.eventBegin(logEvent)

        from pylith.mpi.Communicator import mpi_is_root

        if not interfaces is None:
            cohesiveLabelValue = 100
            for material in problem.materials.components():
                labelValue = material.labelValue
                cohesiveLabelValue = max(cohesiveLabelValue, labelValue+1)
            for interface in interfaces:
                if mpi_is_root():
                    self._info.log("Adjusting topology for fault '%s'." % interface.labelName)
                interface.preinitialize(problem)
                interface.setCohesiveLabelValue(cohesiveLabelValue)
                interface.adjustTopology(mesh)
                cohesiveLabelValue += 1

        self._eventLogger.eventEnd(logEvent)

    def _setupLogging(self):
        """Setup event logging.
        """
        if not "_loggingPrefix" in dir(self):
            self._loggingPrefix = "PL.MeshGenerator."

        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("MeshGenerator")
        logger.initialize()

        events = ["create",
                  "adjTopo"]
        for event in events:
            logger.registerEvent(f"{self._loggingPrefix}{event}")

        self._eventLogger = logger


# End of file
