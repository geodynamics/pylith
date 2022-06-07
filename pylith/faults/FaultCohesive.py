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
from .faults import FaultCohesive as ModuleFaultCohesive


def validateLabel(value):
    """Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError(
            "Label for fault group/nodeset/pset in mesh not specified.")
    return value


def validateDir(value):
    """Validate direction.
    """
    msg = "Direction must be a 3 component vector (list)."
    if not isinstance(value, list):
        raise ValueError(msg)
    if 3 != len(value):
        raise ValueError(msg)
    try:
        nums = list(map(float, value))
    except:
        raise ValueError(msg)
    return nums


class FaultCohesive(Physics, ModuleFaultCohesive):
    """
    Abstract base class for a fault surface implemeted with cohesive cells.
    """

    import pythia.pyre.inventory

    labelName = pythia.pyre.inventory.str("label", default="", validator=validateLabel)
    labelName.meta['tip'] = "Name of label identifier for fault."

    labelValue = pythia.pyre.inventory.int("label_value", default=1)
    labelValue.meta['tip'] = "Value of label identifier for fault."

    edgeName = pythia.pyre.inventory.str("edge", default="")
    edgeName.meta['tip'] = "Name of label identifier for buried fault edges."

    edgeValue = pythia.pyre.inventory.int("edge_value", default=1)
    edgeValue.meta['tip'] = "Value of label identifier for buried fault edges."

    refDir1 = pythia.pyre.inventory.list("ref_dir_1", default=[0.0, 0.0, 1.0], validator=validateDir)
    refDir1.meta['tip'] = "First choice for reference direction to discriminate among tangential directions in 3-D."

    refDir2 = pythia.pyre.inventory.list("ref_dir_2", default=[0.0, 1.0, 0.0], validator=validateDir)
    refDir2.meta['tip'] = "Second choice for reference direction to discriminate among tangential directions in 3-D."

    def __init__(self, name="fault"):
        """Constructor.
        """
        Physics.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Setup fault.
        """
        Physics.preinitialize(self, problem)

        ModuleFaultCohesive.setSurfaceLabelName(self, self.labelName)
        ModuleFaultCohesive.setSurfaceLabelValue(self, self.labelValue)
        ModuleFaultCohesive.setBuriedEdgesLabelName(self, self.edgeName)
        ModuleFaultCohesive.setBuriedEdgesLabelValue(self, self.edgeValue)
        ModuleFaultCohesive.setRefDir1(self, self.refDir1)
        ModuleFaultCohesive.setRefDir2(self, self.refDir2)
        return

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        return

    def _configure(self):
        """Setup members using inventory.
        """
        Physics._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to corresponding C++ object.
        """
        raise NotImplementedError(
            "Please implement _createModuleObj() in derived class.")


# End of file
