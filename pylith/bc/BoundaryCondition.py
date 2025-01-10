# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.problems.Physics import Physics
from .bc import BoundaryCondition as ModuleBoundaryCondition


def validateLabel(value):
    """Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for boundary condition group/nodeset/pset in mesh not specified.")
    return value


class BoundaryCondition(Physics, ModuleBoundaryCondition):
    """
    Abstract base class for boundary conditions.
    """

    import pythia.pyre.inventory

    field = pythia.pyre.inventory.str("field", default="displacement")
    field.meta['tip'] = "Solution subfield associated with boundary condition."

    labelName = pythia.pyre.inventory.str("label", default="", validator=validateLabel)
    labelName.meta['tip'] = "Name of label identifying boundary."

    labelValue = pythia.pyre.inventory.int("label_value", default=1)
    labelValue.meta['tip'] = "Value of label identifying boundary (tag of physical group in Gmsh files)."

    def __init__(self, name="boundarycondition"):
        """Constructor.
        """
        Physics.__init__(self, name, facility="boundary_condition")
        return

    def preinitialize(self, problem):
        """Setup boundary condition.
        """
        Physics.preinitialize(self, problem)

        ModuleBoundaryCondition.setSubfieldName(self, self.field)
        ModuleBoundaryCondition.setLabelName(self, self.labelName)
        ModuleBoundaryCondition.setLabelValue(self, self.labelValue)
        return

    def _configure(self):
        """Setup members using inventory.
        """
        Physics._configure(self)
        return

# End of file
