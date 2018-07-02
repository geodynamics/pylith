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
# @file pylith/meshio/Observer.py
#
# @brief Python abstract base class for observers of subjects (e.g., output).
#
# Factory: observer

from pylith.utils.PetscComponent import PetscComponent
from .feassemble import Observer as ModuleObserver


class Observer(PetscComponent, ModuleObserver):
    """
    Python abstract base class for observers of subjects (e.g., output).

    FACTORY: observer
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="observer"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="observer")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        PetscComponent._configure(self)
        return


# End of file
