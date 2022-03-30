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
        ModuleOutputTriggerTime.setTimeSkip(self, self.timeSkip)

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
