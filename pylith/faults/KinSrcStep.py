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
