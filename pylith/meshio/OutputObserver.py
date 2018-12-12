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
# @file pylith/meshio/OutputObserver.py
#
# @brief Python class for managing output of finite-element information.
#
# Factory: observer

from pylith.utils.PetscComponent import PetscComponent
from .meshio import OutputObserver as ModuleOutputObserver


class OutputObserver(PetscComponent, ModuleOutputObserver):
    """
    Python abstract base class for managing output of finite-element
    information.

    INVENTORY

    Properties
      - None

    Facilities
      - *trigger* Trigger for defining how often output is written.
      - *writer* Writer for output.
      - *field_filter* Filter for output fields.

    FACTORY: observer
    """

    import pyre.inventory

    from .OutputTriggerStep import OutputTriggerStep
    trigger = pyre.inventory.facility("trigger", family="output_trigger", factory=OutputTriggerStep)
    trigger.meta['tip'] = "Trigger defining how often output is written."

    from .DataWriterHDF5 import DataWriterHDF5
    writer = pyre.inventory.facility("writer", factory=DataWriterHDF5, family="data_writer")
    writer.meta['tip'] = "Writer for data."

    from .FieldFilterNone import FieldFilterNone
    fieldFilter = pyre.inventory.facility("field_filter", family="output_field_filter", factory=FieldFilterNone)
    fieldFilter.meta['tip'] = "Filter for output fields."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputmanager"):
        """
        Constructor.
        """
        Observer.__init__(self, name)
        return

    def preinitialize(self, problem):
        """
        Setup output manager.
        """
        self._createModuleObj(observers)
        ModuleOutputObserver.setIdentifier(self, self.aliases[-1])

        self.trigger.preinitialize()
        ModuleOutputObserver.trigger(self, self.trigger)

        self.fieldFilter.preinitialize()
        ModuleOutputObserver.fieldFilter(self, self.fieldFilter)

        self.writer.preinitialize()
        ModuleOutputObserver.writer(self, self.writer)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        Observer._configure(self)
        return

    def _createModuleObj(self, observers):
        """
        Create handle to C++ object.
        """
        raise NotImplementedError("Implement in subclass.")
        return


# End of file
