# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .OutputTrigger import OutputTrigger
from .meshio import OutputTriggerStep as ModuleOutputTriggerStep


class OutputTriggerStep(OutputTrigger, ModuleOutputTriggerStep):
    """
    Define how often output is written in terms of solution steps.

    Implements `OutputTrigger`.
    """
    DOC_CONFIG = {
        "cfg": """
            [output_trigger]
            num_skip = 2
        """
    }

    import pythia.pyre.inventory

    numSkip = pythia.pyre.inventory.int("num_skip", default=0, validator=pythia.pyre.inventory.greaterEqual(0))
    numSkip.meta['tip'] = "Number of solution steps to skip between writes (0 means write every time step)."

    def __init__(self, name="outputtriggerstep"):
        """Constructor.
        """
        OutputTrigger.__init__(self, name)

    def preinitialize(self):
        """Setup output trigger.
        """
        ModuleOutputTriggerStep.__init__(self)
        ModuleOutputTriggerStep.setIdentifier(self, self.aliases[-1])
        ModuleOutputTriggerStep.setNumStepsSkip(self, self.numSkip)

    def _configure(self):
        """Set members based using inventory.
        """
        OutputTrigger._configure(self)

# FACTORIES ////////////////////////////////////////////////////////////


def output_trigger():
    """Factory associated with OutputTriggerStep.
    """
    return OutputTriggerStep()


# End of file
