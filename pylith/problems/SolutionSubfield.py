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
# @file pylith/problems/SolutionSubfield.py
#
# @brief Python object for defining attributes of a subfield within a
# field.
#
# Factory: soln_subfield.

from pylith.topology.Subfield import Subfield


def validateAlias(value):
    """Validate user alias for subfield.
    """
    if 0 == len(value):
        raise ValueError("User-specified alias for subfield not specified.")
    import re
    if re.search(r"\s", value):
        raise ValueError("User-specified alias for subfield cannot contain whitespace.")
    return value


class SolutionSubfield(Subfield):
    """Python object for defining attributes of a subfield within a field.

    FACTORY: soln_subfield
    """

    import pythia.pyre.inventory

    # Override userAlias in derived class with appropriate default.
    userAlias = pythia.pyre.inventory.str("alias", default="", validator=validateAlias)
    userAlias.meta['tip'] = "Name for subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="solution_subfield"):
        """Constructor.
        """
        Subfield.__init__(self, name)

        # Set in derived class initialize().
        self.fieldComponents = None
        self.vectorFieldType = None
        self.scale = None
        self.isFaultOnly = False
        return

    def initialize(self, normalizer, spaceDim):
        """Initialize subfield metadata.
        """
        raise NotImplementedError("Implement in derived class.")

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        from pylith.topology.topology import FieldBase

        Subfield._configure(self)
        return

    def _setComponents(self, spaceDim):
        from pylith.topology.Field import Field
        self.componentNames = []
        if self.vectorFieldType == Field.SCALAR:
            self.componentNames = [self.fieldName]
        elif self.vectorFieldType == Field.VECTOR:
            labels = ["x", "y", "z"]
            self.componentNames = ["{}_{}".format(self.userAlias, label) for label in labels[:spaceDim]]
        else:
            raise NotImplementedError("Not implemented for vector field type %d" % self.vectorFieldType)
        return


# ITEM FACTORIES ///////////////////////////////////////////////////////

def subfieldFactory(name):
    """Factory for subfield items.
    """
    from pythia.pyre.inventory import facility
    return facility(name, family="soln_subfield", factory=SolutionSubfield)


# FACTORIES ////////////////////////////////////////////////////////////

def soln_subfield():
    """Factory associated with SolutionSubfield.
    """
    return SolutionSubfield()


# End of file
