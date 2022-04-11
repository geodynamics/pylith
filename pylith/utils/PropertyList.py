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
