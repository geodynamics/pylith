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

# @file pylith/problems/SolutionSubfield.py
##
# @brief Python object for defining attributes of a subfield within a
# field.
##
# Factory: soln_subfield.

from pylith.utils.PetscComponent import PetscComponent


def validateName(value):
    """
    Validate name of subfield.
    """
    if 0 == len(value):
        raise ValueError("Name of subfield not specified.")
    import re
    if re.search(r"\s", value):
        raise ValueError("Name of subfield cannot contain whitespace.")
    return value


# SolutionSubfield class
class SolutionSubfield(PetscComponent):
    """
    Python object for defining attributes of a subfield within a field.

    Factory: subfield.
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(PetscComponent.Inventory):
        """
        Python object for managing SolutionSubfield facilities and properties.
        """

        # @class Inventory
        # Python object for managing SolutionSubfield facilities and properties.
        ##
        # \b Properties
        # @li \b name Name for subfield.
        # @li \b basis_order Order of basis functions.
        # @li \b quadrature_order Order of numerical quadrature.
        # @li \b is_basis_continuous Is basis continuous?
        # @li \b feSpace Finite-element space [polynomial, point).
        ##
        # \b Facilities
        # @li None

        import pyre.inventory

        # Override in derived class with appropriate default.
        fieldName = pyre.inventory.str("name", default="", validator=validateName)
        fieldName.meta['tip'] = "Name for subfield."

        basisOrder = pyre.inventory.int("basis_order", default=1)
        basisOrder.meta['tip'] = "Order of basis functions."

        quadOrder = pyre.inventory.int("quadrature_order", default=1)
        quadOrder.meta['tip'] = "Order of numerical quadrature."

        isBasisContinuous = pyre.inventory.bool("is_basis_continous", default=True)
        isBasisContinuous.meta['tip'] = "Is basis continuous?"

        feSpaceStr = pyre.inventory.str("finite_element_spave", default="polynomial", validator=pyre.inventory.choice(["polynomial", "point"]))
        feSpaceStr.meta['tip'] = "Finite-element space (polynomial or point). Point space corresponds to delta functions at quadrature points."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="subfield"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="solution_subfield")

        # Set in derived class initialize().
        self.fieldComponents = None
        self.vectorFieldType = None
        self.scale = None
        return

    def initialize(self, normalizer, spaceDim):
        """
        Initialize subfield metadata.
        """
        raise NotImplementedError("Implement in derived class.")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        from pylith.topology.topology import FieldBase

        PetscComponent._configure(self)
        self.fieldName = self.inventory.fieldName
        self.basisOrder = self.inventory.basisOrder
        self.quadOrder = self.inventory.quadOrder
        self.isBasisContinuous = self.inventory.isBasisContinuous
        spaceMapping = {
            "polynomial": FieldBase.POLYNOMIAL_SPACE,
            "point": FieldBase.POINT_SPACE,
        }
        self.feSpace = spaceMapping[self.inventory.feSpaceStr]
        return

    def _setComponents(self, spaceDim):
        from pylith.topology.Field import Field
        self.componentNames = []
        if self.vectorFieldType == Field.SCALAR:
            self.componentNames = self.fieldName
        elif self.vectorFieldType == Field.VECTOR:
            labels = ["x", "y", "z"]
            for iDim in xrange(spaceDim):
                self.componentNames.append("%s_%s" % (self.fieldName, labels[iDim]))
        else:
            raise NotImplementedError("Not implemented for vector field type %d" % self.vectorFieldType)
        return


# ITEM FACTORIES ///////////////////////////////////////////////////////

def subfieldFactory(name):
    """
    Factory for subfield items.
    """
    from pyre.inventory import facility
    return facility(name, family="soln_subfield", factory=SolutionSubfield)


# FACTORIES ////////////////////////////////////////////////////////////

def soln_subfield():
    """
    Factory associated with Subfield.
    """
    return SolutionSubfield()


# End of file
