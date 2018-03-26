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

# @file pylith/bc/Neumann.py
##
# @brief Python object for managing a Neummann (natural) boundary condition.
##
# Factory: boundary_condition

from .BoundaryCondition import BoundaryCondition
from .bc import NeumannTimeDependent as ModuleNeumann
from pylith.feassemble.IntegratorPointwise import IntegratorPointwise

# Neumann class


class Neumann(BoundaryCondition,
                   IntegratorPointwise,
                  ModuleNeumann):
    """
    Python object for managing a Neumann (natural)
    boundary condition.

    Factory: boundary_condition
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(BoundaryCondition.Inventory, IntegratorPointwise.Inventory):
        """
        Python object for managing Neumann facilities and properties.
        """

        # @class Inventory
        # Python object for managing Neumann facilities and properties.
        ##
        # \b Properties
        # @li \b scale Type of scale for nondimenaionlizing Neumann boundary condition (e.g., "pressure" for elasticity").
        ##
        # \b Facilities
        # @li None

        import pyre.inventory

        scaleName = pyre.inventory.str("scale_name", default="pressure", validator=pyre.inventory.choice(["length", "time", "pressure", "density", "velocity"]))
        scaleName.meta['tip'] = "Type of scale for nondimensionalizing Neumann boundary condition ('pressure' for elasticity)."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="neumann"):
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

    def finalize(self):
        """
        Cleanup after running problem.
        """
        print(":TODO: @brad Implement Neumann.finalize() once output manager is added.")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        try:
            BoundaryCondition._configure(self)
            IntegratorPointwise._configure(self)
            ModuleNeumann.scaleName(self, self.inventory.scaleName)
        except ValueError as err:
            aliases = ", ".join(self.aliases)
            raise ValueError("Error while configuring Dirichlet boundary condition "
                             "(%s):\n%s" % (aliases, err.message))
        return

# End of file
