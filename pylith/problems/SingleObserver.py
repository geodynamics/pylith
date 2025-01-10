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


class SingleSolnObserver(PetscComponent):
    """
    Container of solution observers with one observer.
    """

    import pythia.pyre.inventory

    from pylith.meshio.OutputSolnDomain import OutputSolnDomain
    output = pythia.pyre.inventory.facility("domain", family="observer", factory=OutputSolnDomain)
    output.meta['tip'] = "Observer of solution."

    def __init__(self, name="singlesolnobserver"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="singlesolnobserver")


class SinglePhysicsObserver(PetscComponent):
    """
    Container of physics observers with one observer.
    """

    import pythia.pyre.inventory

    from pylith.meshio.OutputPhysics import OutputPhysics
    output = pythia.pyre.inventory.facility("observer", family="observer", factory=OutputPhysics)
    output.meta['tip'] = "Observer of subject."

    def __init__(self, name="singlephysicsobserver"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="singlephysicsobserver")


# End of file
