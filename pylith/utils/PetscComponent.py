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
# @file pylith/utils/PetscComponent.py
#
# @brief Python PetscComponent object for aid in deallocating data
# structures before calling PetscFinalize().

from pyre.components.Component import Component
from .utils import PetscComponent as ModulePetscComponent


class PetscComponent(Component, ModulePetscComponent):
    """
    Python PetscComponent object for aid in deallocating data structures
    before calling PetscFinalize().
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name, facility):
        """
        Constructor.
        """
        Component.__init__(self, name, facility)
        self._createModuleObj()
        return

    def cleanup(self):
        """
        Deallocate data structures.
        """
        for component in self.components():
            if isinstance(component, PetscComponent):
                component.cleanup()

            # Facility arrays are not PetscComponents but have components().
            components = getattr(component, "components", None)
            if callable(components):
                for subcomponent in components():
                    if isinstance(subcomponent, PetscComponent):
                        subcomponent.cleanup()

        self._cleanup()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Create handle to C++ object.
        """
        raise NotImplementedError("Implement in subclass.")

    def _cleanup(self):
        """
        Deallocate locally managed data structures.
        """
        deallocate = getattr(self, "deallocate", None)
        if callable(deallocate):
            deallocate()
        return


# End of file
