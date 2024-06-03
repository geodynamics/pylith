# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .MeshRefiner import MeshRefiner
from .topology import RefineUniform as ModuleRefineUniform


class RefineUniform(MeshRefiner, ModuleRefineUniform):
    """
    Uniform global mesh refinement in parallel.

    Implements `MeshRefiner`.
    """
    DOC_CONFIG = {
        "cfg": """
            # Refine mesh twice to reduce size of cell edges by a factor of 4.
            [pylithapp.mesh_generator.refiner]
            levels = 2
        """
    }

    import pythia.pyre.inventory

    levels = pythia.pyre.inventory.int("levels", default=1, validator=pythia.pyre.inventory.greaterEqual(1))
    levels.meta['tip'] = "Number of refinement levels."

    def __init__(self, name="refineuniform"):
        """Constructor.
        """
        MeshRefiner.__init__(self, name)

    def preinitialize(self):
        """Do minimal initialization."""
        MeshRefiner.preinitialize(self)

        self._createModuleObj()

    def refine(self, mesh):
        """Refine mesh.
        """
        self._setupLogging()
        logEvent = f"{self._loggingPrefix}refine"
        self._eventLogger.eventBegin(logEvent)

        from pylith.mpi.Communicator import mpi_is_root
        if mpi_is_root():
            self._info.log("Refining mesh using uniform refinement.")

        from .Mesh import Mesh
        newMesh = Mesh()
        newMesh.setCoordSys(mesh.getCoordSys())
        ModuleRefineUniform.refine(self, newMesh, mesh, self.levels)
        mesh.cleanup()

        self._eventLogger.eventEnd(logEvent)
        return newMesh

    def _configure(self):
        """Set members based using inventory.
        """
        MeshRefiner._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleRefineUniform.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_refiner():
    """Factory associated with RefineUniform.
    """
    return RefineUniform()


# End of file
