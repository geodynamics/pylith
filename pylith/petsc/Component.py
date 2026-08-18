# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pythia.pyre.components.Component import Component as PyreComponent


class Component(PyreComponent):
    """
    Extension of Pyre Component object for deallocating data structures before calling PetscFinalize().
    """

    def __init__(self, name, facility):
        """Constructor.
        """
        PyreComponent.__init__(self, name, facility)

    def cleanup(self):
        """Deallocate data structures.
        """
        for component in self.components():
            if isinstance(component, Component):
                component.cleanup()

            # Facility arrays are not Components but have components().
            elif hasattr(component, "components"):
                for subcomponent in component.components():
                    if isinstance(subcomponent, Component):
                        subcomponent.cleanup()

        self._cleanup()

    def _cleanup(self):
        """Deallocate locally managed data structures.
        """
        # If module object not yet created, return
        if getattr(self, "this", None) is None:
            return

        deallocate = getattr(self, "deallocate", None)
        if callable(deallocate):
            deallocate()


# End of file
