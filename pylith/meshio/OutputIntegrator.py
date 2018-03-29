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

# @file pyre/meshio/OutputIntegrator.py
##
# @brief Python object for managing output of finite-element
# solution information.
##
# Factory: output_manager

from .OutputManager import OutputManager
from .meshio import OutputIntegrator as ModuleOutputIntegrator

# OutputIntegrator class


class OutputIntegrator(OutputManager, ModuleOutputIntegrator):
    """
    Python object for managing output of finite-element solution
    information.

    @class Inventory
    Python object for managing OutputIntegrator facilities and properties.

    \b Properties
    @li \b vertex_data_fields Names of vertex data fields to output.

    \b Facilities
    @li None

    Factory: output_manager
    """

    # INVENTORY //////////////////////////////////////////////////////////

    class Inventory(OutputManager.Inventory):
        """Python object for managing OutputIntegrator facilities and properties.
        """

        import pyre.inventory

        vertexInfoFields = pyre.inventory.list("vertex_info_fields", default=["all"])
        vertexInfoFields.meta['tip'] = "Names of vertex information fields to output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="OutputIntegrator"):
        """
        Constructor.
        """
        OutputManager.__init__(self, name)
        return

    def preinitialize(self):
        """
        Do
        """
        OutputManager.preinitialize(self)
        ModuleOutputIntegrator.vertexInfoFields(self, self.vertexInfoFields)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        OutputManager._configure(self)
        self.vertexInfoFields = self.inventory.vertexInfoFields
        return

    def _createModuleObj(self):
        """
        Create handle to C++ object.
        """
        ModuleOutputIntegrator.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
    """
    Factory associated with OutputManager.
    """
    return OutputIntegrator()


# End of file
