# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .OutputTrigger import OutputTrigger
from .meshio import OutputTriggerTime as ModuleOutputTriggerTime


class OutputTriggerTime(OutputTrigger, ModuleOutputTriggerTime):
    """
    Define how often output is written in terms of elasped simulation time.

    :::{tip}
    Due to floating point roundoff, it is usually a good idea to use a value that is a fraction of a time step smaller than the desired value.
    :::

    Implements `OutputTrigger`.
    """
    DOC_CONFIG = {
        "cfg": """
            [output_trigger]
            elapsed_time = 0.9999*year
        """
    }

    import pythia.pyre.inventory

    from pythia.pyre.units.time import s
    timeSkip = pythia.pyre.inventory.dimensional("elapsed_time", default=0.0*s)
    timeSkip.meta['tip'] = "Elapsed time between writes."

    def __init__(self, name="outputtriggertime"):
        """Constructor.
        """
        OutputTrigger.__init__(self, name)

    def preinitialize(self):
        """Setup output trigger.
        """
        ModuleOutputTriggerTime.__init__(self)
        ModuleOutputTriggerTime.setIdentifier(self, self.aliases[-1])
        ModuleOutputTriggerTime.setTimeSkip(self, self.timeSkip.value)

    def _configure(self):
        """Set members based using inventory.
        """
        OutputTrigger._configure(self)

# FACTORIES ////////////////////////////////////////////////////////////


def output_trigger():
    """Factory associated with OutputTriggerTime.
    """
    return OutputTriggerTime()


# End of file
