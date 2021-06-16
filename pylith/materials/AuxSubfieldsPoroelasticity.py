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
# @file pylith/materials/AuxSubieldsPoroelasticity.py
#
# @brief Python container for poroelasticity equation subfields.

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsPoroelasticity(PetscComponent):
    """Python container for poroelasticity equation subfields.

    FACTORY: auxiliary_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    porosity = pythia.pyre.inventory.facility("porosity", family="auxiliary_subfield", factory=Subfield)
    porosity.meta['tip'] = "Porosity subfield."

    solidDensity = pythia.pyre.inventory.facility("solid_density", family="auxiliary_subfield", factory=Subfield)
    solidDensity.meta['tip'] = "Solid density subfield."

    fluidDensity = pythia.pyre.inventory.facility("fluid_density", family="auxiliary_subfield", factory=Subfield)
    fluidDensity.meta['tip'] = "Fluid density subfield."

    fluidViscosity = pythia.pyre.inventory.facility("fluid_viscosity", family="auxiliary_subfield", factory=Subfield)
    fluidViscosity.meta['tip'] = "Fluid viscosity subfield."

    bodyForce = pythia.pyre.inventory.facility("body_force", family="auxiliary_subfield", factory=Subfield)
    bodyForce.meta['tip'] = "Body force subfield."

    sourceDensity = pythia.pyre.inventory.facility("source_density", family="auxiliary_subfield", factory=Subfield)
    sourceDensity.meta['tip'] = "Source density subfield."    

    constantPressureSource = pythia.pyre.inventory.facility("constant_pressure_source", family="auxiliary_subfield", factory=Subfield)
    constantPressureSource.meta['tip'] = "Constant pressure source subfield."    

    gravitationalAcceleration = pythia.pyre.inventory.facility(
        "gravitational_acceleration", family="auxiliary_subfield", factory=Subfield)
    gravitationalAcceleration.meta['tip'] = "Gravitational acceleration subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxsubfieldsporoelasticity"):
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
    """Factory associated with AuxSubfieldsPoroelasticity.
    """
    return AuxSubfieldsPoroelasticity()


# End of file
