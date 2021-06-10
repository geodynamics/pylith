# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/materials/SingleProblem.py
#
# @brief Python problems container with one problem.

from pylith.utils.PetscComponent import PetscComponent


class SingleProblem(PetscComponent):
    """Python problems container with one problem.
    """

    import pythia.pyre.inventory

    from .TimeDependent import TimeDependent
    problem = pythia.pyre.inventory.facility("problem", family="problem", factory=TimeDependent)
    problem.meta['tip'] = "Single problem to solve."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="singleproblem"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="problem")
        return


# End of file
