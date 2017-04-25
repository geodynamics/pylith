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

# @file pylith/bc/BoundaryConditionNew.py
##
# @brief Python abstract base class for managing a boundary condition.
##
# This implementation of a boundary condition applies to a single
# boundary of an domain.
##
# Factory: boundary_condition

from pylith.utils.PetscComponent import PetscComponent
from .bc import BoundaryConditionNew as ModuleBoundaryCondition

# Validator for label


def validateLabel(value):
    """
    Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for boundary condition group/nodeset/pset in mesh not specified.")
    return value


# BoundaryConditionNew class
class BoundaryConditionNew(PetscComponent,
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
        Python object for managing BoundaryConditionNew facilities and properties.
        """

        # @class Inventory
        # Python object for managing BoundaryConditionNew facilities and properties.
        ##
        # \b Properties
        # @li \b label Label identifier for boundary.
        ##
        # \b Facilities

        import pyre.inventory

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
        return

    def finalize(self):
        """
        Cleanup after running problem.
        """
        raise NotImplementedError("finalize() not implemented.")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        try:
            PetscComponent._configure(self)
            self.label = self.inventory.label
        except ValueError as err:
            aliases = ", ".join(self.aliases)
            raise ValueError("Error while configuring boundary condition "
                             "(%s):\n%s" % (aliases, err.message))

        return

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError, \
            "Please implement _createModuleObj() in derived class."


# End of file
