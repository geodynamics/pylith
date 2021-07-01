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
# @file pythia.pyre/meshio/OutputSoln.py
#
# @brief Python object for managing output of finite-element
# solution information.
#
# Factory: observer

from .OutputObserver import OutputObserver
from .meshio import OutputSoln as ModuleOutputSoln


class OutputSoln(OutputObserver, ModuleOutputSoln):
    """Python object for managing output of finite-element solution
    information.

    FACTORY: observer
    """

    import pythia.pyre.inventory

    dataFields = pythia.pyre.inventory.list("data_fields", default=["all"])
    dataFields.meta['tip'] = "Names of data fields to output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputsoln"):
        """Constructor.
        """
        OutputObserver.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Do mimimal initialization.
        """
        OutputObserver.preinitialize(self, problem)
        ModuleOutputSoln.setOutputSubfields(self, self.dataFields)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        OutputObserver._configure(self)
        return

# End of file
