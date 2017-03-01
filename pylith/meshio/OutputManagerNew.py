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
from meshio import OutputManagerNew as ModuleOutputManager

# OutputManager class


class OutputManagerNew(PetscComponent, ModuleOutputManager):
    """
    Python abstract base class for managing output of finite-element
    information.

    \b Properties
    @li \b trigger Flag indicating whether to use 'elapsed_time' or 'skip_timesteps' to set how often output is written.
    @li \b elapsed_time Elapsed time between writes.
    @li \b skip_timesteps Number of time steps to skip between writes.

    \b Facilities
    @li \b coordsys Coordinate system for output. NOT IMPLEMENTED
    @li \b vertex_filter Filter for vertex data.
    @li \b cell_filter Filter for cell data.
    """

    # INVENTORY //////////////////////////////////////////////////////////

    import pyre.inventory

    trigger = pyre.inventory.str("trigger", default="num_timesteps", validator=pyre.inventory.choice(["elasped_time", "skip_timesteps"]))
    trigger.meta['tip'] = "Flag indicating whether to use 'elapsed_time' or 'skip_timesteps' to set how often output is written."

    from pyre.units.time import s
    timeSkip = pyre.inventory.dimensional("elapsed_time", default=1.0 * s)
    timeSkip.meta['tip'] = "Elapsed time between writes.""

    numTimeStepsSkip = pyre.inventory.int("skip_timesteps", default=0, validator=pyre.inventory.greaterEqual(0))
    numTimeStepsSkip.meta['tip'] = "Number of time steps to skip between writes."

    from DataWriterHDF5 import DataWriterHDF5
    writer = pyre.inventory.facility("writer", factory=DataWriterHDF5, family="data_writer")
    writer.meta['tip'] = "Writer for data."

    # Not implemented
    #from spatialdata.geocoords.CSCart import CSCart
    #coordsys = pyre.inventory.facility("coordsys", family="coordsys", factory=CSCart)
    #coordsys.meta['tip'] = "Coordinate system for output."

    vertexFilter = pyre.inventory.facility("vertex_filter", family="output_vertex_filter", factory=NullComponent)
    vertexFilter.meta['tip'] = "Filter for vertex data."

    cellFilter = pyre.inventory.facility("cell_filter", family="output_cell_filter", factory=NullComponent)
    cellFilter.meta['tip'] = "Filter for cell data."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputmanager"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="outputmanager")
        self._loggingPrefix = "OutM "
        self._createModuleObj()
        self.coordsys = None  # not implemented
        return

    def preinitialize(self):
        """
        Setup output manager.
        """
        ModuleOutputManager.trigger(self.trigger)
        if self.trigger == ModuleOutputManager.SKIP_TIMESTEPS:
            ModuleOutputManager.numTimeStepsSkip(self.numTimeStepsSkip)
        elif self.trigger == ModuleOutputManager.ELAPSED_TIME:
            ModuleOutputManager.timeSkip(self.timeSkip)
        if self.coordsys:
            ModuleOutputManager.coordsys(self, self.coordsys)
        if self.vertexFilter:
            ModuleOutputManager.vertexFilter(self, self.vertexFilter)
        if self.cellFilter:
            ModuleOutputManager.cellFilter(self, self.cellFilter)
        ModuleOutputManager.writer(self, self.writer)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        PetscComponent._configure(self)
        if self.inventory.trigger == "skip_timesteps":
            self.trigger = ModuleOutputManager.SKIP_TIMESTEPS
            self.numTimeStepsSkip = self.inventory.numTimeStepsSkip
        elif self.inventory.trigger == "elapsed_time":
            self.trigger = ModuleOutputManager.ELAPSED_TIME
            self.timeSkip = self.inventory.timeSkip
        else:
            raise ValueError("Unknown output trigger '%s'." % self.inventory.trigger)

        #self.coordsys = self.inventory.coordsys
        if isinstance(self.inventory.vertexFilter, NullComponent):
            self.vertexFilter = None
        else:
            self.vertexFilter = self.inventory.vertexFilter
        if isinstance(self.inventory.cellFilter, NullComponent):
            self.cellFilter = None
        else:
            self.cellFilter = self.inventory.cellFilter
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
    return OutputManagerMew()


# End of file
