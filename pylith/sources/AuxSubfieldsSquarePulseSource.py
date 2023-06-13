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
# @file pylith/sources/AuxSubfieldsSquarePulseSource.py
#
# @brief Python container for squarepulsesource equation subfields.

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsSquarePulseSource(PetscComponent):
    """Python container for squarepulsesource equation subfields.

    FACTORY: auxiliary_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    volumeFlowRate = pythia.pyre.inventory.facility("volume_flow_rate", family="auxiliary_subfield", factory=Subfield)
    volumeFlowRate.meta['tip'] = "Volume Flow Rate subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxsubfieldssquarepulsesource"):
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
    """Factory associated with AuxSubfieldsSquarePulseSource.
    """
    return AuxSubfieldsSquarePulseSource()


# End of file
