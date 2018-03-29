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

# @file pylith/bc/BoundaryCondition.py
##
# @brief Python abstract base class for managing a boundary condition.
##
# This implementation of a boundary condition applies to a single
# boundary of an domain.
##
# Factory: boundary_condition

from pylith.utils.PetscComponent import PetscComponent
from .bc import BoundaryCondition as ModuleBoundaryCondition

# Validator for label


def validateLabel(value):
    """
    Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for boundary condition group/nodeset/pset in mesh not specified.")
    return value


# BoundaryCondition class
class BoundaryCondition(PetscComponent,
                        ModuleBoundaryCondition):
    """
    Python abstract base class for managing a boundary condition.

    This implementation of a boundary condition applies to a single
    face of an domain.

    Factory: boundary_condition
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(PetscComponent.Inventory):
        """
        Python object for managing BoundaryCondition facilities and properties.
        """

        # @class Inventory
        # Python object for managing BoundaryCondition facilities and properties.
        #
        # \b Properties
        # @li \b label Label identifier for boundary.
        # @li \b field Field in solution to constrain.
        #
        # \b Facilities

        import pyre.inventory

        field = pyre.inventory.str("field", default="displacement")
        field.meta['tip'] = "Solution field associated with boundary condition."

        label = pyre.inventory.str("label", default="", validator=validateLabel)
        label.meta['tip'] = "Label identifier for boundary."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="boundarycondition"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="boundary_condition")
        self._createModuleObj()
        return

    def preinitialize(self, mesh):
        """
        Setup boundary condition.
        """
        ModuleBoundaryCondition.label(self, self.label)
        ModuleBoundaryCondition.field(self, self.field)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        PetscComponent._configure(self)
        self.label = self.inventory.label
        self.field = self.inventory.field
        return

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Please implement _createModuleObj() in derived class.")


# End of file
