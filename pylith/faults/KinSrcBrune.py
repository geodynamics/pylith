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
# @file pylith/faults/KinSrcBrune.py
#
# @brief Python object for a Brune's (1970) slip time function.
#
# Factory: eq_kinematic_src

from .KinSrc import KinSrc
from .faults import KinSrcBrune as ModuleKinSrc


class KinSrcBrune(KinSrc, ModuleKinSrc):
    """Python object for Brune's (1970) far-field slip time function.

    Factory: eq_kinematic_src
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

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
