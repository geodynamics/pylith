# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsFault(PetscComponent):
    """
    Auxiliary subfields associated with a fault.
    """
    DOC_CONFIG = {
        "cfg": """
            # We set the basis order to represent linear variations in the slip subfield.
            [pylithapp.problem.interfaces.fault.auxiliary_fields]
            slip.basis_order = 1
        """
    }

    import pythia.pyre.inventory
    from pylith.topology.Subfield import Subfield

    slip = pythia.pyre.inventory.facility("slip", family="auxiliary_subfield", factory=Subfield)
    slip.meta['tip'] = "Slip subfield."

    def __init__(self, name="auxsubfieldfault"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")

    def _configure(self):
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """Factory associated with AuxSubfieldsFault.
    """
    return AuxSubfieldsFault()


# End of file
