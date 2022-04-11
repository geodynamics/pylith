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

from .KinSrc import KinSrc
from .faults import KinSrcConstRate as ModuleKinSrc


class KinSrcConstRate(KinSrc, ModuleKinSrc):
    """
    Constant slip rate slip time function.

    Implements `KinSrc`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
            origin_time = 100*year

            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.label = Constant slip rate slip time function auxiliary field spatial database
            db_auxiliary_field.values = [initiation_time, slip_rate_left_lateral, slip_rate_opening]
            db_auxiliary_field.data = [0.0*s, -2.0*mm/year, 0.0*mm/year]
            """
    }

    def __init__(self, name="kinsrcconstrate"):
        """Constructor.
        """
        KinSrc.__init__(self, name)
        return

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleKinSrc.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def eq_kinematic_src():
    """Factory associated with KinSrcConstRate.
    """
    return KinSrcConstRate()


# End of file
