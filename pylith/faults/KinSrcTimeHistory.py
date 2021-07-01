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
#
# @file pylith/faults/KinSrcTimeHistory.py
#
# @brief Python object for time history data file slip time function.
#
# Factory: eq_kinematic_src

from .KinSrc import KinSrc
from .faults import KinSrcTimeHistory as ModuleKinSrc


class KinSrcTimeHistory(KinSrc, ModuleKinSrc):
    """Python object for time history data file slip time function.

    Factory: eq_kinematic_src
    """
    import pythia.pyre.inventory

    from spatialdata.spatialdb.TimeHistory import TimeHistory
    dbTimeHistory = pythia.pyre.inventory.facility("time_history", factory=TimeHistory, family="temporal_database")
    dbTimeHistory.meta['tip'] = "Time history with normalized amplitude as a function of time."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="kinsrctimehistory"):
        """Constructor.
        """
        KinSrc.__init__(self, name)
        return

    def preinitialize(self):
        """Do pre-initialization setup.
        """
        KinSrc.preinitialize(self)

        ModuleKinSrc.setTimeHistory(self, self.dbTimeHistory)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleKinSrc.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def eq_kinematic_src():
    """Factory associated with KinSrcTimeHistory.
    """
    return KinSrcTimeHistory()


# End of file
