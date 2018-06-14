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
# @file pylith/feassemble/IntegratorBoundary.py
#
# @brief Python abstract base class for pointwise integrators.

from .IntegratorPointwise import IntegratorPointwise
from .feassemble import IntegratorBoundary as ModuleIntegratorBoundary


def validateLabel(value):
    """
    Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for boundary condition group/nodeset/pset in mesh not specified.")
    return value


def validateDir(value):
    """
    Validate direction.
    """
    msg = "Direction must be a 3 component vector (list)."
    if not isinstance(value, list):
        raise ValueError(msg)
    if 3 != len(value):
        raise ValueError(msg)
    try:
        nums = map(float, value)
    except:
        raise ValueError(msg)
    return nums


class IntegratorBoundary(IntegratorPointwise, ModuleIntegratorBoundary):
    """
    Python abstract base class for integrators over external boundaries.

    INVENTORY

    Properties
      - *label* Label identifying boundary.
      - *field* Field associated with boundary condition.
      - *ref_dir_1* First choice for reference direction to discriminate among tangential directions in 3-D.
      - *ref_dir_2* Second choice for reference direction to discriminate among tangential directions in 3-D.

    Facilities
      - None

    FACTORY: N/A
    """

    import pyre.inventory

    field = pyre.inventory.str("field", default="displacement")
    field.meta['tip'] = "Solution field associated with boundary condition."

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for boundary."

    refDir1 = pyre.inventory.list("ref_dir_1", default=[0.0, 0.0, 1.0], validator=validateDir)
    refDir1.meta['tip'] = "First choice for reference direction to discriminate among tangential directions in 3-D."

    refDir2 = pyre.inventory.list("ref_dir_2", default=[0.0, 1.0, 0.0], validator=validateDir)
    refDir2.meta['tip'] = "Second choice for reference direction to discriminate among tangential directions in 3-D."

# PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="integratorboundary"):
        """
        Constructor.
        """
        IntegratorPointwise.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        IntegratorPointwise.preinitialize(self, mesh)

        ModuleIntegratorBoundary.label(self, self.label)
        ModuleIntegratorBoundary.field(self, self.field)
        ModuleIntegratorBoundary.refDir1(self, self.refDir1)
        ModuleIntegratorBoundary.refDir2(self, self.refDir2)
        return

# PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        IntegratorPointwise._configure(self)
        return


# End of file
