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
# @file pylith/sources/WellboreSource.py
#
# @brief Python object for solving the wellboresource equation.
#
# Factory: source

from .Source import Source
from .sources import WellboreSource as ModuleWellboreSource


class WellboreSource(Source, ModuleWellboreSource):
    """Python source property manager.

    FACTORY: source
    """

    import pythia.pyre.inventory


    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="wellboresource"):
        """Constructor.
        """
        Source.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsWellboreSource import AuxSubfieldsWellboreSource
        self.auxiliarySubfields = AuxSubfieldsWellboreSource("auxiliary_subfields")


    def preinitialize(self, problem):
        """Setup source.
        """
        Source.preinitialize(self, problem)


        return

    def _createModuleObj(self):
        """Create handle to C++ WellboreSource.
        """
        ModuleWellboreSource.__init__(self)
        return


# Factories

def source():
    """Factory associated with WellboreSource.
    """
    return WellboreSource()


# End of file
