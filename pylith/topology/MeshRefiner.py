# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
