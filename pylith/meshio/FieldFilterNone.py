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
# @file pythia.pyre/meshio/FieldFilterNone.py
#
# @brief Python class for null field filter.
#
# Factory: output_field_filter

from .FieldFilter import FieldFilter
from .meshio import FieldFilterNone as ModuleFieldFilterNone


class FieldFilterNone(FieldFilter, ModuleFieldFilterNone):
    """
    Python class for null field filter.

    FACTORY: output_field_filter
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="fieldfilternone"):
        """
        Constructor.
        """
        FieldFilter.__init__(self, name)
        return

    # PRIVATE METHODS /////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleFieldFilterNone.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def output_field_filter():
    """
    Factory associated with FieldFilter.
    """
    return FieldFilterNone()


# End of file
