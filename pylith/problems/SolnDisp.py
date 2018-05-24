#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/problems/SolnDisp.py
##
# @brief Python subfields container with displacement subfield.

from pylith.utils.PetscComponent import PetscComponent


class SolnDisp(PetscComponent):
    """
    Python subfields container with displacement subfield.
    """

    # INVENTORY //////////////////////////////////////////////////////////
    #
    # \b Properties
    # @li None
    #
    # \b Facilities
    # @li \b displacement Displacement subfield.

    import pyre.inventory

    from SubfieldDisplacement import SubfieldDisplacement
    displacement = pyre.inventory.facility("displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="solndisp"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")
        return

    def _configure(self):
        PetscComponent._configure(self)
        return


# End of file
