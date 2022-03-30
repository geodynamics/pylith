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

from .OutputObserver import OutputObserver
from .meshio import OutputSoln as ModuleOutputSoln


class OutputSoln(OutputObserver, ModuleOutputSoln):
    """
    Abstract base class for output of solution subfields.

    Implements `OutputObserver`.
    """

    import pythia.pyre.inventory

    dataFields = pythia.pyre.inventory.list("data_fields", default=["all"])
    dataFields.meta['tip'] = "Names of solution subfields to include in output."

    def __init__(self, name="outputsoln"):
        """Constructor.
        """
        OutputObserver.__init__(self, name)

    def preinitialize(self, problem):
        """Do mimimal initialization.
        """
        OutputObserver.preinitialize(self, problem)
        ModuleOutputSoln.setOutputSubfields(self, self.dataFields)

    def _configure(self):
        """Set members based using inventory.
        """
        OutputObserver._configure(self)

# End of file
