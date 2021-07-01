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
# @file pylith/problems/SolnDisp.py
#
# @brief Python solution container with displacement subfield.

from pylith.utils.PetscComponent import PetscComponent


class SolnDisp(PetscComponent):
    """Python solution container with displacement subfield.
    """

    import pythia.pyre.inventory

    from .SubfieldDisplacement import SubfieldDisplacement
    displacement = pythia.pyre.inventory.facility("displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="solndisp"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")
        return

    def _configure(self):
        PetscComponent._configure(self)
        return

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to insure order is [displacement].
        """
        return [self.displacement]


# Solution uses SolnDisp as the factory, so we do not need override the default
# factory here like we do in the other solution field containers.

# End of file
