# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from .FaultRheology import FaultRheology
from .faults import FrictionStatic as ModuleFrictionStatic


class FrictionStatic(FaultRheology, ModuleFrictionStatic):
    """
    Static friction fault constitutive model.

    Implements `FaultRheology`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.faults.fault.fault_rheology]

            auxiliary_subfields.cohesion.basis_order = 0
            auxiliary_subfields.static_coefficient.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    def __init__(self, name="frictionstatic"):
        """Constructor.
        """
        FaultRheology.__init__(self, name)

    def _defaults(self):
        from .AuxSubfieldsFrictionStatic import AuxSubfieldsFrictionStatic
        self.auxiliarySubfields = AuxSubfieldsFrictionStatic("auxiliary_subfields")

    def preinitialize(self, problem):
        FaultRheology.preinitialize(self, problem)
        ModuleFrictionStatic.useReferenceState(self, self.useReferenceState)

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleFrictionStatic.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def fault_rheology():
    """Factory associated with FrictionStatic.
    """
    return FrictionStatic()


# End of file
