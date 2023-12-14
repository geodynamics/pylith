# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
