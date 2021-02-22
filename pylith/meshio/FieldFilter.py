# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pythia.pyre/meshio/FieldFilter.py
#
# @brief Python abstract base class for filtering vertex fields when
# writing finite-element data.
#
# Factory: output_vertex_filter

from pylith.utils.PetscComponent import PetscComponent


class FieldFilter(PetscComponent):
    """
    Python abstract base class for filtering fields when writing
    finite-element data.

    FACTORY: output_field_filter
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="fieldfilter"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="fieldfilter")
        return

    def preinitialize(self):
        """Do minimal initialization."""
        self._createModuleObj()
        self.setIdentifier(self.aliases[-1])
        return


# End of file
