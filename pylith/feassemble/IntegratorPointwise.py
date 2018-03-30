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
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/feassemble/IntegratorPointwise.py
##
# @brief Python abstract base class for pointwise integrators.

from pylith.utils.PetscComponent import PetscComponent
from .feassemble import IntegratorPointwise as ModuleIntegrator

# IntegratorPointwise class


class IntegratorPointwise(PetscComponent,
                          ModuleIntegrator):
    """
    Python abstract base class for pointwise integrators.

    Factory: material
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(PetscComponent.Inventory):
        """
        Python object for managing IntegratorPointwise facilities and properties.
        """

        # @class Inventory
        # Python object for managing IntegratorPointwise facilities and properties.
        ##
        # \b Properties
        # @li None
        ##
        # \b Facilities
        # @li \b auxiliary_fields Discretization of auxiliary fields associated with material.
        # @li \b db_auxiliary_fields Database for auxiliary fields associated with material.

        import pyre.inventory

        from pylith.topology.AuxSubfield import subfieldFactory
        from pylith.utils.EmptyBin import EmptyBin
        auxSubfields = pyre.inventory.facilityArray("auxiliary_subfields", itemFactory=subfieldFactory, factory=EmptyBin)
        auxSubfields.meta['tip'] = "Discretization of physical properties and state variables."

        from spatialdata.spatialdb.SimpleDB import SimpleDB
        auxFieldDB = pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=SimpleDB)
        auxFieldDB.meta['tip'] = "Database for physical property parameters."

        from pylith.meshio.OutputManager import OutputManager
        outputManager = pyre.inventory.facility("output", family="output_manager", factory=OutputManager)
        outputManager.meta['tip'] = "Output manager."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="integratorpointwise"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="integrator")
        self._createModuleObj()
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        ModuleIntegrator.identifier(self, self.aliases[-1])
        ModuleIntegrator.auxFieldDB(self, self.auxFieldDB)
        ModuleIntegrator.output(self, self.outputManager)

        for subfield in self.auxSubfields.components():
            fieldName = subfield.aliases[-1]
            ModuleIntegrator.auxSubfieldDiscretization(self, fieldName, subfield.basisOrder, subfield.quadOrder, subfield.isBasisContinuous, subfield.feSpace)

        self.outputManager.preinitialize()
        return

# PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        PetscComponent._configure(self)
        self.auxSubfields = self.inventory.auxSubfields
        self.auxFieldDB = self.inventory.auxFieldDB
        self.outputManager = self.inventory.outputManager
        return

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Please implement _createModuleOb() in derived class.")


# End of file
