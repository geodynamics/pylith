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
# @file pylith/meshio/OutputTrigger.py
#
# @brief Python abstract base class for defining how often output is written.
#
# Factory: output_trigger

from pylith.utils.PetscComponent import PetscComponent
from .meshio import OutputTrigger as ModuleOutputTrigger


class OutputTrigger(PetscComponent, ModuleOutputTrigger):
    """Python abstract base class for managing how often output is written.

    FACTORY: output_trigger
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputtrigger"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="outputtrigger")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)
        return


# End of file
