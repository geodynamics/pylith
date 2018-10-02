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
from pylith.problems import Physics as ModulePhysics


class Physics(PetscComponent, ModulePhysics):
    """
    Python abstract base class for objects defining physics.

    INVENTORY

    Properties
      - None

    Facilities
      - *auxiliary_subfields* Discretization of physical properties and state variables.
      - *db_auxiliary_field* Database for physical property parameters.
      - *observers* Observers of integrator (e.g., output).

    FACTORY: material
    """

    import pyre.inventory

    from pylith.topology.AuxSubfield import subfieldFactory
    from pylith.utils.EmptyBin import EmptyBin
    auxiliarySubfields = pyre.inventory.facilityArray(
        "auxiliary_subfields", itemFactory=subfieldFactory, factory=EmptyBin)
    auxiliarySubfields.meta['tip'] = "Discretization of physical properties and state variables."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    auxiliaryFieldDB = pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=SimpleDB)
    auxiliaryFieldDB.meta['tip'] = "Database for physical property parameters."

    observers = pyre.inventory.facilityArray("observers", itemFactory=observerFactory, factory=EmptyBin)
    observers.meta['tip'] = "Observers (e.g., output)."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="physics")
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        self._createModuleObj()

        self.observers.preinitialize(self)

        ModulePhysics.setAuxiliaryFieldDB(self, self.auxiliaryFieldDB)

        for subfield in self.auxiliarySubfields.components():
            fieldName = subfield.aliases[-1]
            ModulePhysics.setAuxiliarySubfieldDiscretization(self, fieldName, subfield.basisOrder, subfield.quadOrder,
                                                             subfield.isBasisContinuous, subfield.feSpace)

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
