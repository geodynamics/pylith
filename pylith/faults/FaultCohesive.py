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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/faults/FaultCohesive.py
#
# @brief Python abstract base class for a fault surface implemented
# with cohesive elements.
#
# Factory: fault

from pylith.feassemble.IntegratorPointwise import IntegratorPointwise
from .faults import FaultCohesive as ModuleFaultCohesive


def validateLabel(value):
    """
    Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for fault group/nodeset/pset in mesh not specified.")
    return value


def validateDir(value):
    """
    Validate direction.
    """
    msg = "Direction must be a 3 component vector (list)."
    if not isinstance(value, list):
        raise ValueError(msg)
    if 3 != len(value):
        raise ValueError(msg)
    try:
        nums = map(float, value)
    except:
        raise ValueError(msg)
    return nums


class FaultCohesive(IntegratorPointwise, ModuleFaultCohesive):
    """
    Python abstract base class for a fault surface implemeted with
    cohesive elements.

    INVENTORY

    Properties
      - *id* Fault identifier
      - *label* Label identifier for fault.
      - *edge* Label identifier for buried fault edges.
      - *ref_dir_1* First choice for reference direction to discriminate among tangential directions in 3-D.
      - *ref_dir_2* Second choice for reference direction to discriminate among tangential directions in 3-D.

    Facilities
      - None

    FACTORY: fault
    """

    # INVENTORY //////////////////////////////////////////////////////////

    import pyre.inventory

    matId = pyre.inventory.int("id", default=100)
    matId.meta['tip'] = "Fault identifier (must be unique across all faults and materials)."

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for fault."

    edge = pyre.inventory.str("edge", default="")
    edge.meta['tip'] = "Label identifier for buried fault edges."

    refDir1 = pyre.inventory.list("ref_dir_1", default=[0.0, 0.0, 1.0], validator=validateDir)
    refDir1.meta['tip'] = "First choice for reference direction to discriminate among tangential directions in 3-D."

    refDir2 = pyre.inventory.list("ref_dir_2", default=[0.0, 1.0, 0.0], validator=validateDir)
    refDir2.meta['tip'] = "Second choice for reference direction to discriminate among tangential directions in 3-D."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="fault"):
        """
        Constructor.
        """
        IntegratorPointwise.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Setup fault.
        """
        ModuleFaultCohesive.id(self, self.matId)
        ModuleFaultCohesive.label(self, self.label)
        ModuleFaultCohesive.edge(self, self.edge)
        ModuleFaultCohesive.refDir1(self, self.refDir1)
        ModuleFaultCohesive.refDir2(self, self.refDir2)
        return

    def verifyConfiguration(self):
        """
        Verify compatibility of configuration.
        """
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        IntegratorPointwise._configure(self)
        return

    def _createModuleObj(self):
        """
        Create handle to corresponding C++ object.
        """
        raise NotImplementedError("Please implement _createModuleObj() in derived class.")


# End of file
