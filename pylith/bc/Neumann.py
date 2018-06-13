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
#
# @brief Python object for managing a Neummann (natural) boundary condition.
#
# Factory: boundary_condition

from .BoundaryCondition import BoundaryCondition
from .bc import NeumannTimeDependent as ModuleNeumann
from pylith.feassemble.IntegratorPointwise import IntegratorPointwise


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


class Neumann(BoundaryCondition, IntegratorPointwise, ModuleNeumann):
    """
    Python object for managing a Neumann (natural)
    boundary condition.

    INVENTORY

    Properties
      - *scale* Type of scale for nondimenaionlizing Neumann boundary condition (e.g., "pressure" for elasticity").

    Facilities
      - *ref_dir_1* First choice for reference direction to discriminate among tangential directions in 3-D.
      - *ref_dir_2* Second choice for reference direction to discriminate among tangential directions in 3-D.

    FACTORY: boundary_condition
    """

    import pyre.inventory

    scaleName = pyre.inventory.str("scale_name", default="pressure",
                                   validator=pyre.inventory.choice(["length", "time", "pressure", "density", "velocity"]))
    scaleName.meta['tip'] = "Type of scale for nondimensionalizing Neumann boundary condition ('pressure' for elasticity)."

    refDir1 = pyre.inventory.list("ref_dir_1", default=[0.0, 0.0, 1.0], validator=validateDir)
    refDir1.meta['tip'] = "First choice for reference direction to discriminate among tangential directions in 3-D."

    refDir2 = pyre.inventory.list("ref_dir_2", default=[0.0, 1.0, 0.0], validator=validateDir)
    refDir2.meta['tip'] = "Second choice for reference direction to discriminate among tangential directions in 3-D."

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

        ModuleNeumann.refDir1(self, self.refDir1)
        ModuleNeumann.refDir2(self, self.refDir2)
        ModuleNeumann.scaleName(self, self.scaleName)
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
