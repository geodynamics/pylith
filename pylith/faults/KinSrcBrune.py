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
from .faults import KinSrcBrune as ModuleKinSrc


class KinSrcBrune(KinSrc, ModuleKinSrc):
    """
    Brune's (1970) far-field slip time function.

    Implements `KinSrc`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
            origin_time = 10*year

            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.label = Brune far-field slip time function auxiliary field spatial database
            db_auxiliary_field.values = [initiation_time, rise_time, final_slip_left_lateral, final_slip_opening]
            db_auxiliary_field.data = [0.0*s, 3.0*s, -2.0*m, 0.0*m]
            """
    }


    def __init__(self, name="kinsrcbrune"):
        """Constructor.
        """
        KinSrc.__init__(self, name)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleKinSrc.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def eq_kinematic_src():
    """Factory associated with KinSrcBrune.
    """
    return KinSrcBrune()


# End of file
