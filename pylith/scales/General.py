# =================================================================================================
# This code is part of SpatialData, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/spatialdata).
#
# Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pythia.pyre.components.Component import Component
from .scales import Scales as ModuleScales


class General(Component, ModuleScales):
    """
    Abstract base class for nondimensionalizing problems.
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="scales"):
        """
        Constructor.
        """
        Component.__init__(self, name, facility="scales")

    def _configure(self):
        Component._configure(self)
        self._createModuleObj()

    def setLengthScale(self, value):
        """
        Set length scale.
        """
        return ModuleScales.setLengthScale(self, value.value)

    def getLengthScale(self):
        """
        Get length scale.
        """
        from pythia.pyre.units.length import meter

        return ModuleScales.getLengthScale(self) * meter

    def setDisplacementScale(self, value):
        """
        Set displacement scale.
        """
        return ModuleScales.setDisplacementScale(self, value.value)

    def getDisplacementScale(self):
        """
        Get displacement scale.
        """
        from pythia.pyre.units.length import meter

        return ModuleScales.getDisplacementScale(self) * meter

    def setRigidityScale(self, value):
        """
        Set pressure scale.
        """
        return ModuleScales.setRigidityScale(self, value.value)

    def getRigidityScale(self):
        """
        Get pressure scale.
        """
        from pythia.pyre.units.pressure import pascal

        return ModuleScales.getRigidityScale(self) * pascal

    def setTimeScale(self, value):
        """
        Get time scale.
        """
        return ModuleScales.setTimeScale(self, value.value)

    def getTimeScale(self):
        """
        Get time scale.
        """
        from pythia.pyre.units.time import second

        return ModuleScales.getTimeScale(self) * second

    def setTemperatureScale(self, value):
        """
        Get temperature scale.
        """
        return ModuleScales.setTemperatureScale(self, value.value)

    def getTemperatureScale(self):
        """
        Get temperature scale.
        """
        from pythia.pyre.units.temperature import kelvin

        return ModuleScales.getTemperatureScale(self) * kelvin

    def nondimensionalize(self, value, scale):
        """
        Make value dimensionless.
        """
        return value / scale

    def dimensionalize(self, value, scale):
        """
        Make value dimensional.
        """
        return value * scale

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Create Python module object.
        """
        ModuleScales.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////


def scales():
    """
    Factory associated with Scales.
    """
    return General()


# End of file
