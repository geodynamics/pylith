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
# @file pylith/problems/Physics.py
#
# @brief Python abstract base class for objects definitng physics, such as behavior of a bulk material,
# boundary condition, interface or constraint.

from pylith.utils.PetscComponent import PetscComponent
from .problems import Physics as ModulePhysics

from pylith.meshio.OutputPhysics import OutputPhysics


# Factories for items in facility arrays

def observerFactory(name):
    """
    Factory for output items.
    """
    from pyre.inventory import facility
    from pylith.meshio.OutputPhysics import OutputPhysics
    return facility(name, family="observer", factory=OutputPhysics)


class Physics(PetscComponent, ModulePhysics):
    """
    Python abstract base class for objects defining physics.

    INVENTORY

    Properties
      - None

    Facilities
      - *db_auxiliary_field* Database for physical property parameters.
      - *observers* Observers of integrator (e.g., output).
    """

    import pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    auxiliaryFieldDB = pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=SimpleDB)
    auxiliaryFieldDB.meta['tip'] = "Database for physical property parameters."

    from pylith.problems.SingleObserver import SinglePhysicsObserver
    observers = pyre.inventory.facilityArray("observers", itemFactory=observerFactory, factory=SinglePhysicsObserver)
    observers.meta['tip'] = "Observers (e.g., output)."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name, facility="physics"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility)

        from pylith.utils.EmptyBin import EmptyBin
        self.auxiliarySubfields = EmptyBin()
        self.derivedSubfields = EmptyBin()
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        self._createModuleObj()
        ModulePhysics.setIdentifier(self, self.aliases[-1])
        ModulePhysics.setAuxiliaryFieldDB(self, self.auxiliaryFieldDB)

        for subfield in self.auxiliarySubfields.components():
            fieldName = subfield.aliases[-1]
            ModulePhysics.setAuxiliarySubfieldDiscretization(self, fieldName, subfield.basisOrder, subfield.quadOrder,
                                                             subfield.dimension, subfield.isBasisContinuous, subfield.feSpace)

        for subfield in self.derivedSubfields.components():
            fieldName = subfield.aliases[-1]
            ModulePhysics.setDerivedSubfieldDiscretization(self, fieldName, subfield.basisOrder, subfield.quadOrder,
                                                           subfield.dimension, subfield.isBasisContinuous, subfield.feSpace)

        for observer in self.observers.components():
            observer.preinitialize(self)
            ModulePhysics.registerObserver(self, observer)
        return

# PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Please implement _createModuleOb() in derived class.")


# End of file
