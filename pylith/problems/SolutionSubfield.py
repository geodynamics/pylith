# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
    """
    Base class for defining attributes of a subfield within a field.
    """

    import pythia.pyre.inventory

    # Set appropriate default in derived class using _defaults().
    userAlias = pythia.pyre.inventory.str("alias", default="", validator=validateAlias)
    userAlias.meta['tip'] = "Name for subfield."

    def __init__(self, name="solution_subfield"):
        """Constructor.
        """
        Subfield.__init__(self, name)

        # Set in derived class initialize().
        self.fieldComponents = None
        self.vectorFieldType = None
        self.scale = None
        self.isFaultOnly = False

    def initialize(self, normalizer, spaceDim):
        """Initialize subfield metadata.
        """
        raise NotImplementedError("Implement in derived class.")

    def _configure(self):
        """Set members based using inventory.
        """
        from pylith.topology.topology import FieldBase
        Subfield._configure(self)

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
