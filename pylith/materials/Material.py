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

from pylith.problems.Physics import Physics
from .materials import Material as ModuleMaterial


def validateLabel(value):
    """Validate descriptive label.
    """
    if 0 == len(value):
        raise ValueError("Descriptive label for material not specified.")
    return value


class Material(Physics, ModuleMaterial):
    """
    Abstract base class for a bulk material.
    """

    import pythia.pyre.inventory

    materialId = pythia.pyre.inventory.int("id", default=0)
    materialId.meta['tip'] = "Material identifier (from mesh generator)."

    label = pythia.pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Descriptive label for material."

    def __init__(self, name="material"):
        """Constructor.
        """
        Physics.__init__(self, name)

    def preinitialize(self, problem):
        """Setup material.
        """
        Physics.preinitialize(self, problem)
        ModuleMaterial.setMaterialId(self, self.materialId)
        ModuleMaterial.setDescriptiveLabel(self, self.label)


# End of file
