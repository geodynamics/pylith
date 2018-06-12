#!/usr/bin/env python
#
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

# @file pyre/meshio/FieldFilterProject.py
#
# @brief Python class for projecting field to another basis.
#
# Factory: output_field_filter

from .FieldFilter import FieldFilter
from .meshio import FieldFilterProject as ModuleFieldFilterProject


class FieldFilterProject(FieldFilter, ModuleFieldFilterProject):
    """
    Python class for projecting field to another basis.

    INVENTORY

    Properties
      - *basis_order* Basis order of projected field.

    Facilities

      - None

    FACTORY: output_field_filter
    """

    import pyre.inventory

    basisOrder = pyre.inventory.int("basis_order", default=1, validator=pyre.inventory.greater(0))
    basisOrder.meta["tip"] = "Basis order of projected field."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="fieldfilterproject"):
        """
        Constructor.
        """
        FieldFilter.__init__(self, name)
        ModuleFieldFilterProject.__init__(self)
        return


    def preinitialize(self):
        """Do minimal initialization."""
        ModuleFieldFilterProject.basisOrder(self, self.basisOrder)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def output_field_filter():
    """
    Factory associated with FieldFilter.
    """
    return FieldFilterProject()


# End of file
