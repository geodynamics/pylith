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
