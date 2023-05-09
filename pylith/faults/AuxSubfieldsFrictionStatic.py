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

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsFrictionStatic(PetscComponent):
    """
    Auxiliary subfields associated with the static friction fault constitutive model.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.faults.fault.fault_rheology.auxiliary_fields]
            cohesion.basis_order = 1
            static_coefficient.basis_order = 1
        """
    }

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    cohesion = pythia.pyre.inventory.facility("cohesion", family="auxiliary_subfield", factory=Subfield)
    cohesion.meta['tip'] = "Cohesion subfield."

    staticCoefficient = pythia.pyre.inventory.facility("static_coefficient", family="auxiliary_subfield", factory=Subfield)
    staticCoefficient.meta['tip'] = "Static coefficient of friction."

    def __init__(self, name="auxsubfieldsfrictionstatic"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")

    def _configure(self):
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """Factory associated with AuxSubfieldsFrictionStatic.
    """
    return AuxSubfieldsFrictionStatic()


# End of file
