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
# @file pylith/materials/AuxSubieldsElasticity.py
#
# @brief Python container for elasticity equation subfields.

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsElasticity(PetscComponent):
    """Python container for elasticity equation subfields.

    FACTORY: auxiliary_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    density = pythia.pyre.inventory.facility("density", family="auxiliary_subfield", factory=Subfield)
    density.meta['tip'] = "Density subfield."

    bodyForce = pythia.pyre.inventory.facility("body_force", family="auxiliary_subfield", factory=Subfield)
    bodyForce.meta['tip'] = "Body force subfield."

    gravitationalAcceleration = pythia.pyre.inventory.facility(
        "gravitational_acceleration", family="auxiliary_subfield", factory=Subfield)
    gravitationalAcceleration.meta['tip'] = "Gravitational acceleration subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxsubfieldselasticity"):
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
    """Factory associated with AuxSubfieldsElasticity.
    """
    return AuxSubfieldsElasticity()


# End of file
