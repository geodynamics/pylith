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
# @file pylith/materials/AuxSubfieldsAbsorbingDampers.py
#
# @brief Python container for absorbing dampers subfields.

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsAbsorbingDampers(PetscComponent):
    """Python container for isotropic, linear elasticity subfields.

    FACTORY: auxiliary_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    density = pythia.pyre.inventory.facility("density", family="auxiliary_subfield", factory=Subfield)
    density.meta['tip'] = "Mass density subfield."

    vs = pythia.pyre.inventory.facility("vs", family="auxiliary_subfield", factory=Subfield)
    vs.meta['tip'] = "Shear (S) wave speed subfield."

    vp = pythia.pyre.inventory.facility("vp", family="auxiliary_subfield", factory=Subfield)
    vp.meta['tip'] = "Dilatational (P) wave speed subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxsubfieldsabsorbingdampers"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        PetscComponent._configure(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """Factory associated with AuxSubfieldsAbsorbingDampers.
    """
    return AuxSubfieldsAbsorbingDampers()


# End of file
