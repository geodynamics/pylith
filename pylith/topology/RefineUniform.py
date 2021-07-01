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
# @file pylith/topology/RefineUniform.py
#
# @brief Python manager for uniform global refinement of mesh in
# parallel.
#
# Factory: mesh_refiner.

from .MeshRefiner import MeshRefiner
from .topology import RefineUniform as ModuleRefineUniform


class RefineUniform(MeshRefiner, ModuleRefineUniform):
    """Python manager for uniform global refinement of mesh in parallel.

    Factory: mesh_refiner
    """

    import pythia.pyre.inventory

    levels = pythia.pyre.inventory.int("levels", default=1, validator=pythia.pyre.inventory.greaterEqual(1))
    levels.meta['tip'] = "Number of refinement levels."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="refineuniform"):
        """Constructor.
        """
        MeshRefiner.__init__(self, name)
        return

    def preinitialize(self):
        """Do minimal initialization."""
        MeshRefiner.preinitialize(self)

        self._createModuleObj()
        return

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
        newMesh.debug(mesh.debug())
        newMesh.setCoordSys(mesh.getCoordSys())
        ModuleRefineUniform.refine(self, newMesh, mesh, self.levels)
        mesh.cleanup()

        self._eventLogger.eventEnd(logEvent)
        return newMesh

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        MeshRefiner._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleRefineUniform.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_refiner():
    """Factory associated with RefineUniform.
    """
    return RefineUniform()


# End of file
