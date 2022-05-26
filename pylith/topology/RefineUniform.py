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
        logEvent = "%srefine" % self._loggingPrefix
        self._eventLogger.eventBegin(logEvent)

        from pylith.mpi.Communicator import petsc_comm_world
        comm = petsc_comm_world()
        if 0 == comm.rank:
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
