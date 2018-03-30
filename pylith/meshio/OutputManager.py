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

# @file pylith/meshio/OutputManager.py
##
# @brief Python abstract base class for managing output of
# finite-element information.
##
# Factory: output_manager

from pylith.utils.PetscComponent import PetscComponent
from pylith.utils.NullComponent import NullComponent
from .meshio import OutputManager as ModuleOutputManager

# OutputManager class


class OutputManager(PetscComponent, ModuleOutputManager):
    """
    Python abstract base class for managing output of finite-element
    information.

    \b Properties
    @li \b trigger Flag indicating whether to use 'elapsed_time' or 'skip_timesteps' to set how often output is written.
    @li \b elapsed_time Elapsed time between writes.
    @li \b skip_timesteps Number of time steps to skip between writes.
    @li \b vertex_info_fields Names of vertex information fields to output.
    @li \b vertex_data_fields Names of vertex data fields to output.
    @li \b cell_info_fields Names of cell information fields to output.
    @li \b cell_data_fields Names of cell data fields to output.

    \b Facilities
    @li \b vertex_filter Filter for vertex data.
    @li \b cell_filter Filter for cell data.
    """

    # INVENTORY //////////////////////////////////////////////////////////

    import pyre.inventory

    trigger = pyre.inventory.str("trigger", default="skip_timesteps", validator=pyre.inventory.choice(["elasped_time", "skip_timesteps"]))
    trigger.meta['tip'] = "Flag indicating whether to use 'elapsed_time' or 'skip_timesteps' to set how often output is written."

    from pyre.units.time import s
    timeSkip = pyre.inventory.dimensional("elapsed_time", default=1.0 * s)
    timeSkip.meta['tip'] = "Elapsed time between writes."

    numTimeStepsSkip = pyre.inventory.int("skip_timesteps", default=0, validator=pyre.inventory.greaterEqual(0))
    numTimeStepsSkip.meta['tip'] = "Number of time steps to skip between writes."

    from DataWriterHDF5 import DataWriterHDF5
    writer = pyre.inventory.facility("writer", factory=DataWriterHDF5, family="data_writer")
    writer.meta['tip'] = "Writer for data."

    vertexFilter = pyre.inventory.facility("vertex_filter", family="output_vertex_filter", factory=NullComponent)
    vertexFilter.meta['tip'] = "Filter for vertex data."

    cellFilter = pyre.inventory.facility("cell_filter", family="output_cell_filter", factory=NullComponent)
    cellFilter.meta['tip'] = "Filter for cell data."

    vertexInfoFields = pyre.inventory.list("vertex_info_fields", default=["all"])
    vertexInfoFields.meta['tip'] = "Names of vertex information fields to output."

    vertexDataFields = pyre.inventory.list("vertex_data_fields", default=["all"])
    vertexDataFields.meta['tip'] = "Names of vertex data fields to output."

    cellInfoFields = pyre.inventory.list("cell_info_fields", default=["all"])
    cellInfoFields.meta['tip'] = "Names of cell information fields to output."

    cellDataFields = pyre.inventory.list("cell_data_fields", default=["all"])
    cellDataFields.meta['tip'] = "Names of cell data fields to output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputmanager"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="outputmanager")
        self._loggingPrefix = "OutM "
        self._createModuleObj()
        return

    def preinitialize(self):
        """
        Setup output manager.
        """
        ModuleOutputManager.identifier(self, self.aliases[-1])
        if self.inventory.trigger == "skip_timesteps":
            ModuleOutputManager.trigger(self, ModuleOutputManager.SKIP_TIMESTEPS)
            ModuleOutputManager.numTimeStepsSkip(self, self.numTimeStepsSkip)
        elif self.trigger == "elapsed_time":
            ModuleOutputManager.trigger(self, ModuleOutputManager.ELAPSED_TIME)
            ModuleOutputManager.timeSkip(self, self.timeSkip)
        else:
            raise ValueError("Unknown output trigger '%s'." % self.inventory.trigger)

        if not isinstance(self.inventory.vertexFilter, NullComponent):
            ModuleOutputManager.vertexFilter(self, self.vertexFilter)
        if not isinstance(self.inventory.cellFilter, NullComponent):
            ModuleOutputManager.cellFilter(self, self.cellFilter)

        ModuleOutputManager.vertexInfoFields(self, self.inventory.vertexInfoFields)
        ModuleOutputManager.vertexDataFields(self, self.inventory.vertexDataFields)
        ModuleOutputManager.cellInfoFields(self, self.inventory.cellInfoFields)
        ModuleOutputManager.cellDataFields(self, self.inventory.cellDataFields)

        self.writer.preinitialize()
        ModuleOutputManager.writer(self, self.writer)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        PetscComponent._configure(self)
        self.trigger = self.inventory.trigger
        self.numTimeStepsSkip = self.inventory.numTimeStepsSkip
        self.timeSkip = self.inventory.timeSkip

        self.coordsys = NullComponent()
        self.writer = self.inventory.writer
        return

    def _createModuleObj(self):
        """
        Create handle to C++ object.
        """
        ModuleOutputManager.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
    """
    Factory associated with OutputManager.
    """
    return OutputManager()


# End of file
