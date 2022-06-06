# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

from pylith.problems.Physics import Physics
from .materials import Material as ModuleMaterial


def validateDescription(value):
    """Validate description.
    """
    if 0 == len(value):
        raise ValueError("Description for material not specified.")
    return value


class Material(Physics, ModuleMaterial):
    """
    Abstract base class for a bulk material.
    """

    import pythia.pyre.inventory

    description = pythia.pyre.inventory.str("description", default="", validator=validateDescription)
    description.meta['tip'] = "Descriptive label for material."

    labelName = pythia.pyre.inventory.str("label", default="material-id", validator=pythia.pyre.inventory.choice(["material-id"]))
    labelName.meta['tip'] = "Name of label for material. Currently only 'material-id' is allowed."

    labelValue = pythia.pyre.inventory.int("label_value", default=1)
    labelValue.meta["tip"] = "Value of label for material."

    def __init__(self, name="material"):
        """Constructor.
        """
        Physics.__init__(self, name)

    def preinitialize(self, problem):
        """Setup material.
        """
        Physics.preinitialize(self, problem)
        ModuleMaterial.setDescription(self, self.description)
        ModuleMaterial.setLabelName(self, self.labelName)
        ModuleMaterial.setLabelValue(self, self.labelValue)


# End of file
