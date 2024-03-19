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
# @file pylith/sources/TimeHistoryWavelet.py
#
# @brief Python source time functiof for a ricker wavelet.
#
# Factory: pointforce_sourcetimefunction

from .SourceTimeFunctionMomentTensorForce import SourceTimeFunctionMomentTensorForce
from .sources import TimeHistoryWavelet as ModuleTimeHistoryWavelet


class TimeHistoryWavelet(SourceTimeFunctionMomentTensorForce, ModuleTimeHistoryWavelet):
    """Python source time function for time history source.

    FACTORY: pointforce_sourcetimefunction
    """

    import pythia.pyre.inventory

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="timehistorywavelet"):
        """Constructor.
        """
        SourceTimeFunctionMomentTensorForce.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsSourceTime import AuxSubfieldsSourceTime
        self.auxiliarySubfields = AuxSubfieldsSourceTime("auxiliary_subfields")

    def preinitialize(self, problem):
        SourceTimeFunctionMomentTensorForce.preinitialize(self, problem)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleTimeHistoryWavelet.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def momenttensorforce_sourcetimefunction():
    """Factory associated with TimeHistoryWavelet.
    """
    return TimeHistoryWavelet()


# End of file
