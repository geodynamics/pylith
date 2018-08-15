# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/meshio/Observers.py
#
# @brief Python abstract base class subjects of observers.
#
# Factory: N/A

from pylith.utils.PetscComponent import PetscComponent
from .feassemble import Observers as ModuleObservers


class Observers(PetscComponent, ModuleObservers):
    """
    Python abstract base class for subjects of observers.

    FACTORY: N/A
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name, facility):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility)
        return

    def preinitialize(self):
        """Do minimal initialization. Mostly setup low-level C++ object."""
        ModuleObservers.identifier(self, self.aliases[-1])
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        PetscComponent._configure(self)
        return


# End of file
