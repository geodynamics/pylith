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
#
# @file pylith/bc/BoundaryCondition.py
#
# @brief Python abstract base class for managing a boundary condition.
#
# This implementation of a boundary condition applies to a single
# boundary of an domain.
#
# Factory: boundary_condition

from pylith.problems.Physics import Physics
from .bc import BoundaryCondition as ModuleBoundaryCondition


def validateLabel(value):
    """Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for boundary condition group/nodeset/pset in mesh not specified.")
    return value


class BoundaryCondition(Physics,
                        ModuleBoundaryCondition):
    """
    Python abstract base class for managing a boundary condition.

    This implementation of a boundary condition applies to a single
    boundary of an domain.

    FACTORY: boundary_condition
    """

    import pythia.pyre.inventory

    field = pythia.pyre.inventory.str("field", default="displacement")
    field.meta['tip'] = "Solution subfield associated with boundary condition."

    label = pythia.pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for boundary."

    def __init__(self, name="boundarycondition"):
        """Constructor.
        """
        Physics.__init__(self, name, facility="boundary_condition")
        return

    def preinitialize(self, problem):
        """Setup boundary condition.
        """
        Physics.preinitialize(self, problem)

        ModuleBoundaryCondition.setMarkerLabel(self, self.label)
        ModuleBoundaryCondition.setSubfieldName(self, self.field)
        return

    def _configure(self):
        """Setup members using inventory.
        """
        Physics._configure(self)
        return

# End of file
