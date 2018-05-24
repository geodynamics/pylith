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
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/bc/AbsorbingDampers.py
##
# @brief Python object for managing an absorbing boundary condition.
##
# Factory: boundary_condition

from .BoundaryCondition import BoundaryCondition
from pylith.feassemble.IntegratorPointwise import IntegratorPointwise


class AbsorbingDampers(BoundaryCondition,
                       IntegratorPointwise):
    """
    Python object for managing an absorbing boundary condition.

    Factory: boundary_condition
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="absorbingdampers"):
        """
        Constructor.
        """
        BoundaryCondition.__init__(self, name)
        IntegratorPointwise.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        BoundaryCondition.preinitialize(self, mesh)
        IntegratorPointwise.preinitialize(self, mesh)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        BoundaryCondition._configure(self)
        IntegratorPointwise._configure(self)
        return

# End of file
