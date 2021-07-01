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
# @file pylith/faults/KinSrcConstRate.py
#
# @brief Python object for a constant slip rate slip time function.
#
# Factory: eq_kinematic_src

from .KinSrc import KinSrc
from .faults import KinSrcConstRate as ModuleKinSrc


class KinSrcConstRate(KinSrc, ModuleKinSrc):
    """Python object for a constant slip rate slip time function.

    Factory: eq_kinematic_src
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="kinsrcconstrate"):
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
    """Factory associated with KinSrcConstRate.
    """
    return KinSrcConstRate()


# End of file
