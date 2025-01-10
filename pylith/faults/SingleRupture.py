# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.utils.PetscComponent import PetscComponent


class SingleRupture(PetscComponent):
    """
    Kinematic slip source container with one source.

    :::{seealso}
    See [`FaultCohesiveKin` Component](FaultCohesiveKin.md).
    :::
    """

    import pythia.pyre.inventory

    from .KinSrcStep import KinSrcStep
    rupture = pythia.pyre.inventory.facility("rupture", family="eq_kinematic_src", factory=KinSrcStep)
    rupture.meta['tip'] = "Kinematic slip source (for example, an earthquake rupture) for fault."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="singlerupture"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="rupture")
        return


# End of file
