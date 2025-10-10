# =================================================================================================
# This code is part of SpatialData, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/spatialdata).
#
# Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from .General import General


class QuasistaticPoroelasticity(General):
    """
    Convenience object for nondimensionalizing quasi-static poroelasticity problems.

    Implements `General`.
    """

    DOC_CONFIG = {
        "cfg": """
            [normalizer]
            length_scale = 100.0*km
            displacement_scale = 50.0*km
            shear_modulus = 25.0*GPa
            viscosity = 0.001*Pa*s
            permeability = 1.0e-12*m**2
            """,
    }

    import pythia.pyre.inventory

    from pythia.pyre.units.pressure import pascal, GPa
    from pythia.pyre.units.length import meter, km
    from pythia.pyre.units.time import year, second

    lengthScale = pythia.pyre.inventory.dimensional("length_scale", default=100.0 * km)
    lengthScale.validator = pythia.pyre.inventory.greater(0.0 * meter)
    lengthScale.meta["tip"] = (
        "Length scale in boundary value problem (size of feature controlling displacement, fault)."
    )

    displacementScale = pythia.pyre.inventory.dimensional(
        "displacement_scale", default=1.0 * meter
    )
    displacementScale.validator = pythia.pyre.inventory.greater(0.0 * meter)
    displacementScale.meta["tip"] = (
        "Nominal displacement scale in boundary value problem."
    )

    shearModulus = pythia.pyre.inventory.dimensional(
        "shear_modulus", default=10.0 * GPa
    )
    shearModulus.validator = pythia.pyre.inventory.greater(0.0 * pascal)
    shearModulus.meta["tip"] = "Nominal shear modulus in boundary value problem."

    viscosity = pythia.pyre.inventory.dimensional(
        "viscosity", default=0.001 * pascal * second
    )
    viscosity.validator = pythia.pyre.inventory.greater(0.0 * pascal)
    viscosity.meta["tip"] = "Nominal fluid viscosity in boundary value problem."

    permeability = pythia.pyre.inventory.dimensional(
        "permeability", default=1.0e-13 * meter**2
    )
    permeability.validator = pythia.pyre.inventory.greater(0.0 * pascal)
    permeability.meta["tip"] = "Nominal permeability in boundary value problem."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="quasistaticporoelasticity"):
        """
        Constructor.
        """
        General.__init__(self, name)

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        from .ElasticityScales import ElasticityScales

        General._configure(self)

        ElasticityScales.setQuasistaticPoroelasticity(
            self,
            lengthScale=self.inventory.lengthScale,
            permeability=self.inventory.permeability,
            viscosity=self.inventory.viscosity,
            rigidity=self.inventory.shearModulus,
        )

        self.setDisplacementScale(self.inventory.displacementScale)


# FACTORIES ////////////////////////////////////////////////////////////


def scales():
    """
    Factory associated with QuasistaticPoroelasticity.
    """
    return QuasistaticPoroelasticity()


# End of file
