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
# @file pylith/meshio/OutputTriggerTime.py
#
# @brief Python class for defining how often output is written in terms of elapsed time.
#
# Factory: output_manager

from .OutputTrigger import OutputTrigger
from .meshio import OutputTriggerTime as ModuleOutputTriggerTime


class OutputTriggerTime(OutputTrigger, ModuleOutputTriggerTime):
    """Python class for defining how often output is writtern in terms of elaspsed time.
    """

    import pythia.pyre.inventory

    from pythia.pyre.units.time import s
    timeSkip = pythia.pyre.inventory.dimensional("elapsed_time", default=0.0*s)
    timeSkip.meta['tip'] = "Elapsed time between writes."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputtriggertime"):
        """Constructor.
        """
        OutputTrigger.__init__(self, name)
        return

    def preinitialize(self):
        """Setup output trigger.
        """
        ModuleOutputTriggerTime.__init__(self)
        ModuleOutputTriggerTime.setIdentifier(self, self.aliases[-1])
        ModuleOutputTriggerTime.setTimeSkip(self, self.timeSkip)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        OutputTrigger._configure(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def output_trigger():
    """Factory associated with OutputTriggerTime.
    """
    return OutputTriggerTime()


# End of file
