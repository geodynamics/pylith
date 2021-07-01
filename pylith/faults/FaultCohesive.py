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
# @file pylith/faults/FaultCohesive.py
#
# @brief Python abstract base class for a fault surface implemented
# with cohesive elements.
#
# Factory: fault

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
    """Python abstract base class for a fault surface implemeted with
    cohesive elements.

    FACTORY: fault
    """

    import pythia.pyre.inventory

    matId = pythia.pyre.inventory.int("id", default=100)
    matId.meta['tip'] = "Fault identifier (must be unique across all faults and materials)."

    label = pythia.pyre.inventory.str(
        "label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for fault."

    edge = pythia.pyre.inventory.str("edge", default="")
    edge.meta['tip'] = "Label identifier for buried fault edges."

    refDir1 = pythia.pyre.inventory.list(
        "ref_dir_1", default=[0.0, 0.0, 1.0], validator=validateDir)
    refDir1.meta['tip'] = "First choice for reference direction to discriminate among tangential directions in 3-D."

    refDir2 = pythia.pyre.inventory.list(
        "ref_dir_2", default=[0.0, 1.0, 0.0], validator=validateDir)
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

        ModuleFaultCohesive.setInterfaceId(self, self.matId)
        ModuleFaultCohesive.setSurfaceMarkerLabel(self, self.label)
        ModuleFaultCohesive.setBuriedEdgesMarkerLabel(self, self.edge)
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
