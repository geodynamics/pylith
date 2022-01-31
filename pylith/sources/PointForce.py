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
# @file pylith/sources/PointForce.py
#
# @brief Python object for solving the pointforce equation.
#
# Factory: source

from pylith.sources.RickerFunction import RickerFunction
from .Source import Source
from .sources import PointForce as ModulePointForce


class PointForce(Source, ModulePointForce):
    """Python source property manager.

    FACTORY: source
    """

    import pythia.pyre.inventory

    source_time_function = pythia.pyre.inventory.facility("source_time_function", family="pointforce_sourcetimefunction", factory=RickerFunction)
    source_time_function.meta['tip'] = "Source time function for point force."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="pointforce"):
        """Constructor.
        """
        Source.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsPointForce import AuxSubfieldsPointForce
        self.auxiliarySubfields = AuxSubfieldsPointForce("auxiliary_subfields")


    def preinitialize(self, problem):
        """Setup source.
        """
        Source.preinitialize(self, problem)


        return

    def _createModuleObj(self):
        """Create handle to C++ PointForce.
        """
        ModulePointForce.__init__(self)
        return


# Factories

def source():
    """Factory associated with PointForce.
    """
    return PointForce()


# End of file
