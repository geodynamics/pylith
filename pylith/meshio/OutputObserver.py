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
# @file pylith/meshio/OutputObserver.py
#
# @brief Python class for managing output of finite-element information.
#
# Factory: observer

from pylith.utils.PetscComponent import PetscComponent
from .meshio import OutputObserver as ModuleOutputObserver


class OutputObserver(PetscComponent, ModuleOutputObserver):
    """Python abstract base class for managing output of finite-element
    information.

    FACTORY: observer
    """

    import pythia.pyre.inventory

    from .OutputTriggerStep import OutputTriggerStep
    trigger = pythia.pyre.inventory.facility(
        "trigger", family="output_trigger", factory=OutputTriggerStep)
    trigger.meta['tip'] = "Trigger defining how often output is written."

    from .DataWriterHDF5 import DataWriterHDF5
    writer = pythia.pyre.inventory.facility(
        "writer", factory=DataWriterHDF5, family="data_writer")
    writer.meta['tip'] = "Writer for data."

    outputBasisOrder = pythia.pyre.inventory.int("output_basis_order", default=1, validator=pythia.pyre.inventory.choice((0,1)))
    outputBasisOrder.meta['tip'] = "Basis order for output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputobserver"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="outputobserver")
        return

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

        self.writer.preinitialize()
        ModuleOutputObserver.setWriter(self, self.writer)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        PetscComponent._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        raise NotImplementedError("Implement in subclass.")
        return


# End of file
