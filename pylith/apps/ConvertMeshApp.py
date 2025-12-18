# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from .PetscApplication import PetscApplication

from pylith.problems.TimeDependent import TimeDependent
from pylith.scales.General import General


class ConvertMeshApp(PetscApplication):
    """Application for pre-initializing mesh (reordering, converting formats)."""

    import pythia.pyre.inventory

    from pylith.initializers.Convert import Initializer

    meshInitializer = pythia.pyre.inventory.facility(
        "mesh_initializer", family="mesh_initializer", factory=Initializer
    )
    meshInitializer.meta["tip"] = "Mesh initializer."

    def __init__(self, name="convertmeshapp"):
        """Constructor."""
        PetscApplication.__init__(self, name)

    def main(self, *args, **kwds):
        """Run the application."""
        scales = General()
        scales._configure()

        problem = TimeDependent()
        problem._createModuleObj()
        problem.setScales(scales)
        self.meshInitializer.preinitialize()
        self.meshInitializer.runPhases(problem)

    def _configure(self):
        """Setup members using inventory."""
        PetscApplication._configure(self)


# End of file
