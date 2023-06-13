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
# @file pylith/sources/MomentTensorForce.py
#
# @brief Python object for solving the momenttensorforce equation.
#
# Factory: source

from pylith.sources.RickerWavelet import RickerWavelet
from .Source import Source
from .sources import MomentTensorForce as ModuleMomentTensorForce


class MomentTensorForce(Source, ModuleMomentTensorForce):
    """Python source property manager.

    FACTORY: source
    """

    import pythia.pyre.inventory

    source_time_function = pythia.pyre.inventory.facility("source_time_function", family="momenttensorforce_sourcetimefunction", factory=RickerWavelet)
    source_time_function.meta['tip'] = "Source time function for momenttensor force."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="momenttensorforce"):
        """Constructor.
        """
        Source.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsMomentTensorForce import AuxSubfieldsMomentTensorForce
        self.auxiliarySubfields = AuxSubfieldsMomentTensorForce("auxiliary_subfields")


    def preinitialize(self, problem):
        """Setup source.
        """
        self.source_time_function.preinitialize(problem)
        Source.preinitialize(self, problem)

        self.source_time_function.addAuxiliarySubfields(self, problem)

        return

    def _createModuleObj(self):
        """Create handle to C++ MomentTensorForce.
        """
        ModuleMomentTensorForce.__init__(self)
        ModuleMomentTensorForce.setSourceTimeFunction(self, self.source_time_function)  # Material sets auxiliary db in source_time_function.
        return


# Factories

def source():
    """Factory associated with MomentTensorForce.
    """
    return MomentTensorForce()


# End of file
