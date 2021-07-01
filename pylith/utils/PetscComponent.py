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
# @file pylith/utils/PetscComponent.py
#
# @brief Python PetscComponent object for aid in deallocating data
# structures before calling PetscFinalize().

from pythia.pyre.components.Component import Component


class PetscComponent(Component):
    """Python PetscComponent object for aid in deallocating data structures
    before calling PetscFinalize().
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name, facility):
        """Constructor.
        """
        Component.__init__(self, name, facility)
        return

    def cleanup(self):
        """Deallocate data structures.
        """
        for component in self.components():
            if isinstance(component, PetscComponent):
                component.cleanup()

            # Facility arrays are not PetscComponents but have components().
            elif hasattr(component, "components"):
                for subcomponent in component.components():
                    if isinstance(subcomponent, PetscComponent):
                        subcomponent.cleanup()

        self._cleanup()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _cleanup(self):
        """Deallocate locally managed data structures.
        """
        # If module object not yet created, return
        if getattr(self, "this", None) is None:
            return

        deallocate = getattr(self, "deallocate", None)
        if callable(deallocate):
            deallocate()
        return


# End of file
