# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pylith.utils.PetscComponent import PetscComponent
from .problems import Physics as ModulePhysics

from pylith.meshio.OutputPhysics import OutputPhysics
from pylith.utils.NullComponent import NullComponent


def observerFactory(name):
    """Factory for output items."""
    from pythia.pyre.inventory import facility
    from pylith.meshio.OutputPhysics import OutputPhysics

    return facility(name, family="observer", factory=OutputPhysics)


class Physics(PetscComponent, ModulePhysics):
    """
    Abstract base class for objects defining physics.
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import subfieldFactory
    from pylith.utils.EmptyBin import EmptyBin

    auxiliarySubfields = pythia.pyre.inventory.facilityArray(
        "auxiliary_subfields", itemFactory=subfieldFactory, factory=EmptyBin
    )
    auxiliarySubfields.meta["tip"] = (
        "Discretization information for auxiliary subfields."
    )

    derivedSubfields = pythia.pyre.inventory.facilityArray(
        "derived_subfields", itemFactory=subfieldFactory, factory=EmptyBin
    )
    derivedSubfields.meta["tip"] = "Discretization of derived subfields."

    from spatialdata.spatialdb.SimpleDB import SimpleDB

    auxiliaryFieldDB = pythia.pyre.inventory.facility(
        "db_auxiliary_field", family="spatial_database", factory=SimpleDB
    )
    auxiliaryFieldDB.meta["tip"] = "Database for physical property parameters."

    from pylith.problems.SingleObserver import SinglePhysicsObserver

    observers = pythia.pyre.inventory.facilityArray(
        "observers", itemFactory=observerFactory, factory=SinglePhysicsObserver
    )
    observers.meta["tip"] = "Observers (e.g., output)."

    def __init__(self, name="physics", facility="physics"):
        """Constructor."""
        PetscComponent.__init__(self, name, facility)

    def preinitialize(self, problem):
        """Do pre-initialization setup."""
        self._createModuleObj()
        identifier = self.aliases[-1]
        ModulePhysics.setIdentifier(self, identifier)
        ModulePhysics.setScales(self, problem.scales)

        if not isinstance(self.auxiliaryFieldDB, NullComponent):
            ModulePhysics.setAuxiliaryFieldDB(self, self.auxiliaryFieldDB)

        for subfield in self.auxiliarySubfields.components():
            fieldName = subfield.aliases[-1]
            descriptor = subfield.getTraitDescriptor("quadrature_order")
            if (
                hasattr(descriptor.locator, "source")
                and descriptor.locator.source == "default"
            ):
                quadOrder = problem.defaults.quadOrder
            else:
                quadOrder = subfield.quadOrder
            ModulePhysics.setAuxiliarySubfieldDiscretization(
                self,
                fieldName,
                subfield.basisOrder,
                quadOrder,
                subfield.dimension,
                subfield.cellBasis,
                subfield.feSpace,
                subfield.isBasisContinuous,
            )

        for subfield in self.derivedSubfields.components():
            fieldName = subfield.aliases[-1]
            descriptor = subfield.getTraitDescriptor("quadrature_order")
            if (
                hasattr(descriptor.locator, "source")
                and descriptor.locator.source == "default"
            ):
                quadOrder = problem.defaults.quadOrder
            else:
                quadOrder = subfield.quadOrder
            ModulePhysics.setDerivedSubfieldDiscretization(
                self,
                fieldName,
                subfield.basisOrder,
                quadOrder,
                subfield.dimension,
                subfield.cellBasis,
                subfield.feSpace,
                subfield.isBasisContinuous,
            )

        for observer in self.observers.components():
            observer.preinitialize(problem, identifier)
            ModulePhysics.registerObserver(self, observer)

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object."""
        raise NotImplementedError(
            "Please implement _createModuleOb() in derived class."
        )


# End of file
