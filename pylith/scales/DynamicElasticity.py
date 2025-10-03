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


class DynamicElasticity(General):
    """
    Convenience object for nondimensionalizing dynamic elasticity problems.

    Implements `General`.
    """

    DOC_CONFIG = {
        "cfg": """
            [normalizer]
            length_scale = 100.0*km
            displacement_scale = 50.0*km
            shear_modulus = 25.0*GPa
            velocity_scale = 3.0*km/s
            """,
    }

    import pythia.pyre.inventory

    from pythia.pyre.units.pressure import pascal, GPa
    from pythia.pyre.units.length import meter, km
    from pythia.pyre.units.time import second

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

    shearWaveSpeed = pythia.pyre.inventory.dimensional(
        "shear_wave_speed", default=3.0 * km / second
    )
    shearWaveSpeed.validator = pythia.pyre.inventory.greater(0.0 * second)
    shearWaveSpeed.meta["tip"] = (
        "Time scale of boundary value problem (for example, viscoelastic relaxation time)."
    )

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="dynamicelasticity"):
        """
        Constructor.
        """
        Scales.__init__(self, name)

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        Scales._configure(self)

        displacement = self.inventory.displacementScale
        length = self.inventory.lengthScale
        shearModulus = self.inventory.shearModulus
        shearWaveSpeed = self.inventory.shearWaveSpeed

        rigidityScale = displacement / length * shearModulus
        timeScale = length / shearWaveSpeed

        self.setLengthScale(length)
        self.setDisplacementScale(displacement)
        self.setRigidityScale(rigidityScale)
        self.setTimeScale(timeScale)


# FACTORIES ////////////////////////////////////////////////////////////


def scales():
    """
    Factory associated with DynamicElasticity.
    """
    return DynamicElasticity()


# End of file
