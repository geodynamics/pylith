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


class AuxSubfieldsElasticity(PetscComponent):
    """
    Auxiliary subfields associated with the elasticity equation.

    Setting the parameters for a subfield does not turn on its use.
    The [`Elasticity` Component](Elasticity.md) has flags for including or excluding terms in the elasticity equation.
    """
    DOC_CONFIG = {
        "cfg": """
            # We set the basis order to represent linear variations in the density and body 
            # force subfields and a uniform gravitational acceleration subfield.
            [pylithapp.problem.materials.mat_elastic.auxiliary_fields]
            density.basis_order = 1
            body_force.basis_order = 1
            gravitational_acceleration.basis_order = 0
        """
    }

    import pythia.pyre.inventory
    from pylith.topology.Subfield import Subfield

    density = pythia.pyre.inventory.facility("density", family="auxiliary_subfield", factory=Subfield)
    density.meta['tip'] = "Density subfield."

    bodyForce = pythia.pyre.inventory.facility("body_force", family="auxiliary_subfield", factory=Subfield)
    bodyForce.meta['tip'] = "Body force subfield."

    gravitationalAcceleration = pythia.pyre.inventory.facility("gravitational_acceleration", family="auxiliary_subfield", factory=Subfield)
    gravitationalAcceleration.meta['tip'] = "Gravitational acceleration subfield."

    def __init__(self, name="auxsubfieldselasticity"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")

    def _configure(self):
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """Factory associated with AuxSubfieldsElasticity.
    """
    return AuxSubfieldsElasticity()


# End of file
