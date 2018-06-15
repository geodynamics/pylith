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
# @file pylith/bc/AbsorbingDampers.py
#
# @brief Python object for managing absorbing dampers condition.
#
# Factory: boundary_condition

from pylith.bc.IntegratorBoundary import IntegratorBoundary
from .bc import AbsorbingDampers as ModuleAbsorbingDampers


class AbsorbingDampers(IntegratorBoundary, ModuleAbsorbingDampers):
    """
    Python object for managing absorbing dampers condition.

    FACTORY: boundary_condition
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="absorbingdampers"):
        """
        Constructor.
        """
        IntegratorBoundary.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        IntegratorBoundary.preinitialize(self, mesh)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        IntegratorBoundary._configure(self)
        return

    def _createModuleObj(self):
        """
        Create handle to corresponding C++ object.
        """
        ModuleAbsorbingDampers.__init__(self)
        return

# End of file
