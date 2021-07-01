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
#
# @file pylith/topology/MeshRefiner.py
#
# @brief Python manager for refining mesh in parallel.
#
# Factory: mesh_refiner.

from pylith.utils.PetscComponent import PetscComponent


class MeshRefiner(PetscComponent):
    """Python abstract base class for refining mesh in parallel.

    Factory: mesh_refiner
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="refiner"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="refiner")
        return

    def preinitialize(self):
        """Do minimal initialization."""
        return

    def refine(self, mesh):
        """Refine mesh.
        """
        return mesh

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)
        return

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
        return


# End of file
