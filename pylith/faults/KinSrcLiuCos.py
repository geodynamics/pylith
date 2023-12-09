# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file pylith/faults/KinSrcLiuCos.py
#
# @brief Python object for a Liu () cosine-sine slip time function.
#
# Factory: eq_kinematic_src

from .KinSrc import KinSrc
from .faults import KinSrcLiuCos as ModuleKinSrc


class KinSrcLiuCos(KinSrc, ModuleKinSrc):
    """
    {cite:t}`Liu:etal:2006` cosine-sine slip time function.

    Implements `KinSrc`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
            origin_time = 10*year

            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Liu cosine-sine slip time function auxiliary field spatial database
            db_auxiliary_field.values = [initiation_time, rise_time, final_slip_left_lateral, final_slip_opening]
            db_auxiliary_field.data = [0.0*s, 3.0*s, -2.0*m, 0.0*m]
            """
    }


    def __init__(self, name="kinsrcliucos"):
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
    """Factory associated with KinSrcLiuCos.
    """
    return KinSrcLiuCos()


# End of file
