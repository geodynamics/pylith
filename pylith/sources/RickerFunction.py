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
# @file pylith/sources/RickerFunction.py
#
# @brief Python source time functiof for a ricker wavelet.
#
# Factory: pointforce_sourcetimefunction

from .SourceTimeFunctionPointForce import SourceTimeFunctionPointForce
from .sources import RickerFunction as ModuleRickerFunction


class RickerFunction(SourceTimeFunctionPointForce, ModuleRickerFunction):
    """Python source time function for ricker source.

    FACTORY: pointforce_sourcetimefunction
    """

    import pythia.pyre.inventory

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="rickerfunction"):
        """Constructor.
        """
        SourceTimeFunctionPointForce.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsRickerFunction import AuxSubfieldsRickerFunction
        self.auxiliarySubfields = AuxSubfieldsRickerFunction("auxiliary_subfields")

    def preinitialize(self, problem):
        SourceTimeFunctionPointForce.preinitialize(self, problem)


        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleRickerFunction.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def pointforce_sourcetimefunction():
    """Factory associated with RickerFunction.
    """
    return RickerFunction()


# End of file
