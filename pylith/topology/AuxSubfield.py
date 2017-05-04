#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/topology/AuxSubfield.py
##
# @brief Python object for defining attributes of a subfield within a
# field.
##
# Factory: subfield.

from pyre.components.Component import Component

# AuxSubfield class


class AuxSubfield(Component):
    """
    Python object for defining discretization of an auxiliary subfield.

    Factory: auxiliary_subfield.
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(Component.Inventory):
        """
        Python object for managing AuxSubfield facilities and properties.
        """

        # @class Inventory
        # Python object for managing AuxSubfield facilities and properties.
        ##
        # \b Properties
        # @li \b basis_order Order of basis functions.
        # @li \b quadrature_order Order of numerical quadrature.
        # @li \b basis_continuous Is basis continuous?
        # @li \b feSpace Finite-element space [polynomial, point).
        ##
        # \b Facilities
        # @li None

        import pyre.inventory

        basisOrder = pyre.inventory.int("basis_order", default=1)
        basisOrder.meta['tip'] = "Order of basis functions."

        quadOrder = pyre.inventory.int("quadrature_order", default=1)
        quadOrder.meta['tip'] = "Order of numerical quadrature."

        isBasisContinuous = pyre.inventory.bool("is_basis_continous", default=True)
        isBasisContinuous.meta['tip'] = "Is basis continuous?"

        feSpaceStr = pyre.inventory.str("finite_element_spave", default="polynomial", validator=pyre.inventory.choice(["polynomial", "point"]))
        feSpaceStr.meta['tip'] = "Finite-element space (polynomial or point). Point space corresponds to delta functions at quadrature points."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxsubfield"):
        """
        Constructor.
        """
        Component.__init__(self, name, facility="auxsubfield")

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        from pylith.topology.topology import FieldBase

        Component._configure(self)
        self.basisOrder = self.inventory.basisOrder
        self.quadOrder = self.inventory.quadOrder
        self.isBasisContinuous = self.inventory.isBasisContinuous
        spaceMapping = {
            "polynomial": FieldBase.POLYNOMIAL_SPACE,
            "point": FieldBase.POINT_SPACE,
        }
        self.feSpace = spaceMapping[self.inventory.feSpaceStr]
        return


# ITEM FACTORIES ///////////////////////////////////////////////////////

def subfieldFactory(name):
    """
    Factory for subfield items.
    """
    from pyre.inventory import facility
    from pylith.topology.AuxSubfield import AuxSubfield
    return facility(name, family="auxiliary_subfield", factory=AuxSubfield)


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfield():
    """
    Factory associated with AuxSubfield.
    """
    return AuxSubfield()


# End of file
