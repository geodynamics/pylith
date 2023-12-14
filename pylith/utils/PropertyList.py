# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pythia.pyre.components.Component import Component


class PropertyList(Component):
    """
    Container for holding a list of items (one per line).
    """

    def __init__(self, name="propertylist"):
        """Constructor.
        """
        Component.__init__(self, name, facility="property_list")
        self.items = []
        return

    def updateConfiguration(self, registry):
        """Update Pyre configuration.
        """
        self.items = [
            (name, descriptor) for name, descriptor in registry.properties.items()
        ]
        return []


# FACTORIES ////////////////////////////////////////////////////////////

def property_list():
    return PropertyList()
