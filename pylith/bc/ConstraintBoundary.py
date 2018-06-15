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
# @file pylith/bc/ConstraintBoundary.py
#
# @brief Python abstract base class for pointwise integrators.

from pylith.feassemble.ConstraintPointwise import ConstraintPointwise
from .bc import ConstraintBoundary as ModuleConstraintBoundary


def validateLabel(value):
    """
    Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for boundary condition group/nodeset/pset in mesh not specified.")
    return value


class ConstraintBoundary(ConstraintPointwise, ModuleConstraintBoundary):
    """
    Python abstract base class for integrators over external boundaries.

    INVENTORY

    Properties
      - *label* Label identifying boundary.
      - *field* Field associated with boundary condition.

    Facilities
      - None

    FACTORY: N/A
    """

    import pyre.inventory

    field = pyre.inventory.str("field", default="displacement")
    field.meta['tip'] = "Solution field associated with boundary condition."

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for boundary."

# PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="constraintboundary"):
        """
        Constructor.
        """
        ConstraintPointwise.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        ConstraintPointwise.preinitialize(self, mesh)

        ModuleConstraintBoundary.label(self, self.label)
        ModuleConstraintBoundary.field(self, self.field)
        return

# PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        ConstraintPointwise._configure(self)
        return


# End of file
