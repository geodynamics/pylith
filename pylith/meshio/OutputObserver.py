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
from .meshio import OutputObserver as ModuleOutputObserver


class OutputObserver(PetscComponent, ModuleOutputObserver):
    """
    Abstract base class for managing output of solution information.
    """

    import pythia.pyre.inventory

    from .OutputTriggerStep import OutputTriggerStep
    trigger = pythia.pyre.inventory.facility("trigger", family="output_trigger", factory=OutputTriggerStep)
    trigger.meta['tip'] = "Trigger defining how often output is written."

    from .DataWriterHDF5 import DataWriterHDF5
    writer = pythia.pyre.inventory.facility("writer", factory=DataWriterHDF5, family="data_writer")
    writer.meta['tip'] = "Writer for data."

    outputBasisOrder = pythia.pyre.inventory.int("output_basis_order", default=1, validator=pythia.pyre.inventory.choice([0,1]))
    outputBasisOrder.meta['tip'] = "Basis order for output."

    refineLevels = pythia.pyre.inventory.int("refine_levels", default=0, validator=pythia.pyre.inventory.greaterEqual(0))
    refineLevels.meta['tip'] = "Number of mesh refinement levels for output."

    def __init__(self, name="outputobserver"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="outputobserver")

    def preinitialize(self, problem):
        """Setup output manager.
        """
        self._createModuleObj()
        ModuleOutputObserver.setIdentifier(self, self.aliases[-1])

        self.trigger.preinitialize()
        ModuleOutputObserver.setTrigger(self, self.trigger)

        descriptor = self.getTraitDescriptor("output_basis_order")
        if hasattr(descriptor.locator, "source") and descriptor.locator.source == "default":
            outputBasisOrder = problem.defaults.outputBasisOrder
        else:
            outputBasisOrder = self.outputBasisOrder
        ModuleOutputObserver.setOutputBasisOrder(self, outputBasisOrder)
        ModuleOutputObserver.setRefineLevels(self, self.refineLevels)

        self.writer.preinitialize()
        ModuleOutputObserver.setWriter(self, self.writer)

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        raise NotImplementedError("Implement in subclass.")


# End of file
