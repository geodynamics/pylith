# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

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
