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
from .faults import KinSrcStep as ModuleKinSrc


class KinSrcStep(KinSrc, ModuleKinSrc):
    """
    Step slip time function.

    Implements `KinSrc`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
            origin_time = 10*year

            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Step slip time function auxiliary field spatial database
            db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
            db_auxiliary_field.data = [0.0*s, -2.0*m, 0.0*m]
            """
    }

    def __init__(self, name="kinsrcstep"):
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
    """Factory associated with KinSrcStep.
    """
    return KinSrcStep()


# End of file
