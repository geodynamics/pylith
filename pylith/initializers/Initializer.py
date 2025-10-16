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

from .initializers import Initializer as ModuleInitializer


def phaseFactory(name):
    """Factory for output items."""
    from pythia.pyre.inventory import facility
    from .InitializePhase import InitializePhase

    return facility(name, family="initialize_phase", factory=InitializePhase)


class Initializer(PetscComponent, ModuleInitializer):
    """
    Manager for reading and setting up a finite-element mesh.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.mesh_initializer]
            phases = pylith.initializers.Serial
            
            # Equivalent manual construction of Serial phases
            phases = [read_mesh, reorder_mesh, distribute_mesh, insert_interfaces, refine_mesh]
            read_mesh = pylith.initializers.MeshReader
            reorder_mesh = pylith.initializers.MeshReordering
            distribute_mesh = pylith.initializers.MeshDistributor
            insert_interfaces = pylith.initializers.InsertInterfaces
            refine_mesh = pylith.initializers.MeshRefiner
        """
    }

    import pythia.pyre.inventory
    from .Serial import Serial

    phases = pythia.pyre.inventory.facilityArray(
        "phases", itemFactory=phaseFactory, factory=Serial
    )
    phases.meta["tip"] = "Phases for mesh initialization."

    def __init__(self, name="mesh_initializer"):
        """Constructor."""
        PetscComponent.__init__(self, name, facility="mesh_initializer")

    def preinitialize(self):
        """Read and setup a finite-element mesh."""
        for phase in self.phases.components():
            phase.preinitialize()
        ModuleInitializer.__init__(self)
        ModuleInitializer.setPhases(self, self.phases.components())

    def _configure(self):
        """Set members based using inventory."""
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////
def mesh_initializer():
    """Factory associated with Initializer."""
    return Initializer()


# End of file
