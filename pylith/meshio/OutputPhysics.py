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
# @file pyre/meshio/OutputPhysics.py
#
# @brief Python object for managing output over points with constrained degrees of freedom.
#
# Factory: output_manager

from .OutputObserver import OutputObserver
from .meshio import OutputPhysics as ModuleOutputPhysics


class OutputPhysics(OutputObserver, ModuleOutputPhysics):
    """
    Python object for managing output over points with constrained degrees of freedom.

    INVENTORY

    Properties
      - *info_fields* Names of info fields to output.
      - *data_fields* Names of data fields to output.

    Facilities
      - None

    Factory: observer
    """

    import pyre.inventory

    infoFields = pyre.inventory.list("info_fields", default=["all"])
    infoFields.meta['tip'] = "Names of info fields to output."

    dataFields = pyre.inventory.list("data_fields", default=["all"])
    dataFields.meta['tip'] = "Names of data fields to output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputphysics"):
        """
        Constructor.
        """
        OutputObserver.__init__(self, name)
        return

    def preinitialize(self, problem):
        """
        Do mimimal initialization.
        """
        OutputObserver.preinitialize(self, problem)
        ModuleOutputSoln.setInfoFieldNames(self, self.infoFields)
        ModuleOutputSoln.setDataFieldNames(self, self.dataFields)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self, integrator):
        """
        Create handle to C++ object.
        """
        ModuleOutputPhysics.__init__(self, integrator)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """
    Factory associated with OutputObserver.
    """
    return OutputPhysics()


# End of file
