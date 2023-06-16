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
