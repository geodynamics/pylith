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
# @file pylith/sources/AuxSubieldsWellboreSource.py
#
# @brief Python container for wellboresource equation subfields.

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsWellboreSource(PetscComponent):
    """Python container for wellboresource equation subfields.

    FACTORY: auxiliary_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    fluidDensity = pythia.pyre.inventory.facility("fluid_density", family="auxiliary_subfield", factory=Subfield)
    fluidDensity.meta['tip'] = "Fluid density subfield."

    fluidViscosity = pythia.pyre.inventory.facility("fluid_viscosity", family="auxiliary_subfield", factory=Subfield)
    fluidViscosity.meta['tip'] = "Fluid viscosity subfield."

    isotropicPermeability = pythia.pyre.inventory.facility(
        "isotropic_permeability", family="auxiliary_subfield", factory=Subfield)
    isotropicPermeability.meta['tip'] = "Isotropic permeability subfield."

    wellboreRadius = pythia.pyre.inventory.facility(
        "wellbore_radius", family="auxiliary_subfield", factory=Subfield)
    wellboreRadius.meta['tip'] = "Wellbore radius subfield."

    wellboreLength = pythia.pyre.inventory.facility(
        "wellbore_length", family="auxiliary_subfield", factory=Subfield)
    wellboreLength.meta['tip'] = "Wellbore length subfield."    

    wellboreCharacter = pythia.pyre.inventory.facility(
        "wellbore_character", family="auxiliary_subfield", factory=Subfield)
    wellboreCharacter.meta['tip'] = "Wellbore character subfield."    

    wellborePressure = pythia.pyre.inventory.facility(
        "wellbore_pressure", family="auxiliary_subfield", factory=Subfield)
    wellborePressure.meta['tip'] = "Wellbore pressure subfield."   

    elementDimensions = pythia.pyre.inventory.facility(
        "element_dimensions", family="auxiliary_subfield", factory=Subfield)
    elementDimensions.meta['tip'] = "Element dimension subfield."    

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxsubfieldswellboresource"):
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
    """Factory associated with AuxSubfieldsWellboreSource.
    """
    return AuxSubfieldsWellboreSource()


# End of file
