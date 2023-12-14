# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .PetscComponent import PetscComponent


class NullComponent(PetscComponent):
    """
    Empty Pyre component.
    """

    def __init__(self):
        """Constructor.
        """
        PetscComponent.__init__(self, name="nullcomponent", facility="nullcomponent")

    def _cleanup(self):
        """Deallocate locally managed data structures.
        """
        return


# End of file
