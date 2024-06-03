# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.bc.BoundaryCondition import BoundaryCondition
from .bc import AbsorbingDampers as ModuleAbsorbingDampers


class AbsorbingDampers(BoundaryCondition, ModuleAbsorbingDampers):
    """
    Absorbing dampers boundary condition.

    Implements `BoundaryCondition`.
    """
    DOC_CONFIG = {
        "cfg": """
            [bc]
            label = boundary_xpos
            field = velocity

            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Material properties for absorbing boundary
            db_auxiliary_field.values = [density, vs, vp]
            db_auxiliary_field.data = [2500*kg/m**3, 1.0*km/s, 1.732*km/s]

            auxiliary_subfields.density.basis_order = 0
            auxiliary_subfields.vp.basis_order = 0
            auxiliary_subfields.vs.basis_order = 0            
            """,
    }

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="absorbingdampers"):
        """Constructor.
        """
        BoundaryCondition.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsAbsorbingDampers import AuxSubfieldsAbsorbingDampers
        self.auxiliarySubfields = AuxSubfieldsAbsorbingDampers("auxiliary_subfields")

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        BoundaryCondition.preinitialize(self, problem)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        BoundaryCondition._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to corresponding C++ object.
        """
        ModuleAbsorbingDampers.__init__(self)
        return


# Factories

def boundary_condition():
    """Factory associated with AbsorbingDampers.
    """
    return AbsorbingDampers()


# End of file
