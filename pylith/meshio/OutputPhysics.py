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
# @file pythia.pyre/meshio/OutputPhysics.py
#
# @brief Python object for managing output over points with constrained degrees of freedom.
#
# Factory: output_manager

from .OutputObserver import OutputObserver
from .meshio import OutputPhysics as ModuleOutputPhysics


class OutputPhysics(OutputObserver, ModuleOutputPhysics):
    """Python object for managing output over points with constrained degrees of freedom.

    Factory: observer
    """

    import pythia.pyre.inventory

    infoFields = pythia.pyre.inventory.list("info_fields", default=["all"])
    infoFields.meta['tip'] = "Names of info fields to output."

    dataFields = pythia.pyre.inventory.list("data_fields", default=["all"])
    dataFields.meta['tip'] = "Names of data fields to output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputphysics"):
        """Constructor.
        """
        OutputObserver.__init__(self, name)
        return

    def preinitialize(self, problem, identifier):
        """Do mimimal initialization.
        """
        OutputObserver.preinitialize(self, problem)
        ModuleOutputPhysics.setInfoFields(self, self.infoFields)
        ModuleOutputPhysics.setDataFields(self, self.dataFields)

        self.writer.setFilename(problem.defaults.outputDir, problem.defaults.simName, identifier)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleOutputPhysics.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """Factory associated with OutputObserver.
    """
    return OutputPhysics()


# End of file
