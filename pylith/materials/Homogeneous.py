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
