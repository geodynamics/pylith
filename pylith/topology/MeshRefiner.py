# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from pylith.utils.PetscComponent import PetscComponent


class MeshRefiner(PetscComponent):
    """
    Abstract base class for refining a mesh in parallel.
    """

    def __init__(self, name="refiner"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="refiner")

    def preinitialize(self):
        """Do minimal initialization."""
        return

    def refine(self, mesh):
        """Refine mesh.
        """
        return mesh

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)

    def _setupLogging(self):
        """Setup event logging.
        """
        self._loggingPrefix = "Refin "
        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("FE Refinement")
        logger.initialize()
        events = ["refine"]
        for event in events:
            logger.registerEvent("%s%s" % (self._loggingPrefix, event))

        self._eventLogger = logger


# End of file
