# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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
            db_auxiliary_field.description = Constant slip rate slip time function auxiliary field spatial database
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
