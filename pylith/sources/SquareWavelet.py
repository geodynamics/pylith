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
# @file pylith/sources/SquareWavelet.py
#
# @brief Python source time functiof for a square wavelet.
#
# Factory: pointforce_sourcetimefunction

from .SourceTimeFunctionMomentTensorForce import SourceTimeFunctionMomentTensorForce
from .sources import SquareWavelet as ModuleSquareWavelet


class SquareWavelet(SourceTimeFunctionMomentTensorForce, ModuleSquareWavelet):
    """Python source time function for square source.

    FACTORY: pointforce_sourcetimefunction
    """

    import pythia.pyre.inventory

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="squarewavelet"):
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
        ModuleSquareWavelet.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def momenttensorforce_sourcetimefunction():
    """Factory associated with SquareWavelet.
    """
    return SquareWavelet()


# End of file
