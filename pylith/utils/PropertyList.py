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
# @file pylith/utils/PropertyList.py
#
# @brief Python class for holding a list of Pyre properties.
#
# The properties are given one per line or one per command line argument.

from pythia.pyre.components.Component import Component


class PropertyList(Component):
    """Python object for holding a list of items (one per line).

    FACTORY: property_list
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
