# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/feassemble/ConstraintPointwise.py
#
# @brief Python abstract base class for constraints on operator
# actions with finite-elements.

from pylith.feassemble.Observers import Observers
from .feassemble import ConstraintPointwise as ModuleConstraint

import numpy


def validateDOF(value):
    """
    Validate list of constrained degrees of freedom.
    """
    try:
        size = len(value)
        num = map(int, value)
        for v in num:
            if v < 0:
                raise ValueError
    except:
        raise ValueError, \
            "'constrained_dof' must be a zero based list of indices of degrees of " \
            "freedom at a vertex."
    return num


def observerFactory(name):
    """
    Factory for output items.
    """
    from pyre.inventory import facility
    from pylith.meshio.OutputConstraint import OutputConstraint
    return facility(name, family="observer", factory=OutputConstraint)


class ConstraintPointwise(Observers, ModuleConstraint):
    """
    Python abstract base class for constraints on operator
    actions with finite-elements.

    INVENTORY

    Properties
      - *constrained_dof* Constrained degrees of freedom.

    Facilities
      - *auxiliary_subfields* Discretization of constraint parameter auxiliary subfields.
      - *db_auxiliary_field* Spatial database for constrain parameters.
      - *observers* Observers of constraint (e.g., output).
    """

    import pyre.inventory

    constrainedDOF = pyre.inventory.list("constrained_dof", default=[], validator=validateDOF)
    constrainedDOF.meta['tip'] = "Constrained degrees of freedom (0=1st DOF, 1=2nd DOF, etc)."

    from pylith.topology.AuxSubfield import subfieldFactory
    from pylith.utils.EmptyBin import EmptyBin
    auxSubfields = pyre.inventory.facilityArray("auxiliary_subfields", itemFactory=subfieldFactory, factory=EmptyBin)
    auxSubfields.meta['tip'] = "Discretization of constraint parameters."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    auxFieldDB = pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=SimpleDB)
    auxFieldDB.meta['tip'] = "Database for constraint parameters."

    from pylith.feassemble.SingleObserver import SingleConstraintObserver
    observers = pyre.inventory.facilityArray("observers", itemFactory=observerFactory, factory=SingleConstraintObserver)
    observers.meta['tip'] = "Observers (e.g., output) for constraint."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="constraintpointwise"):
        """
        Constructor.
        """
        Observers.__init__(self, name, facility="constraint")
        return

    def preinitialize(self, mesh):
        """
        Setup constraint.
        """
        self._createModuleObj()
        Observers.preinitialize(self)

        ModuleConstraint.constrainedDOF(self, numpy.array(self.constrainedDOF, dtype=numpy.int32))
        ModuleConstraint.auxFieldDB(self, self.auxFieldDB)

        for subfield in self.auxSubfields.components():
            fieldName = subfield.aliases[-1]
            ModuleConstraint.auxSubfieldDiscretization(
                self, fieldName, subfield.basisOrder, subfield.quadOrder, subfield.isBasisContinuous, subfield.feSpace)

        for observer in self.observers.components():
            observer.preinitialize(self)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        Observers._configure(self)
        return

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Please implement _createModuleOb() in derived class.")


# End of file
