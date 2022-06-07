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
from .meshio import OutputTrigger as ModuleOutputTrigger


class OutputTrigger(PetscComponent, ModuleOutputTrigger):
    """
    Abstract base class for managing how often output is written.
    """

    def __init__(self, name="outputtrigger"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="outputtrigger")

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)


# End of file
