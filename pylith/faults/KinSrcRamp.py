# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .KinSrc import KinSrc
from .faults import KinSrcRamp as ModuleKinSrc


class KinSrcRamp(KinSrc, ModuleKinSrc):
    """
    Linear ramp slip time function.

    Implements `KinSrc`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
            origin_time = 10*year

            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Ramp slip time function auxiliary field spatial database
            db_auxiliary_field.values = [initiation_time, rise_time, final_slip_left_lateral, final_slip_opening]
            db_auxiliary_field.data = [0.0*s, 3.0*s, -2.0*m, 0.0*m]
            """
    }


    def __init__(self, name="kinsrcramp"):
        """Constructor.
        """
        KinSrc.__init__(self, name)
        return

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleKinSrc.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def eq_kinematic_src():
    """Factory associated with KinSrcRamp.
    """
    return KinSrcRamp()


# End of file
