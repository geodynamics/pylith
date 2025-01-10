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


class Homogeneous(PetscComponent):
    """
    Materials container with one material.

    :::{seealso}
    See [`Problems` Component](../problems/Problem.md).
    :::
    """

    import pythia.pyre.inventory

    from .IsotropicLinearElasticity import IsotropicLinearElasticity
    material = pythia.pyre.inventory.facility("material", family="material", factory=IsotropicLinearElasticity)
    material.meta['tip'] = "Material in problem."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="homogeneous"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="material")


# End of file
