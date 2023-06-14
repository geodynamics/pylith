# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Vuffalo

# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/sources/AuxSubfieldsPointForce.py
#
# @brief Python container for pointforce equation subfields.

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsPointForce(PetscComponent):
    """Python container for pointforce equation subfields.

    FACTORY: auxiliary_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    pointForce = pythia.pyre.inventory.facility("point_force", family="auxiliary_subfield", factory=Subfield)
    pointForce.meta['tip'] = "Point Force subfield."

    TimeDelay = pythia.pyre.inventory.facility("time_delay", family="auxiliary_subfield", factory=Subfield)
    TimeDelay.meta['tip'] = "time delay subfield."


    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxsubfieldspointforce"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        PetscComponent._configure(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """Factory associated with AuxSubfieldsPointForce.
    """
    return AuxSubfieldsPointForce()


# End of file
