#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file tests/topology/test_meshmem.py

## @brief Python application for testing memory
## allocation/deallocation in Mesh.

from pyre.applications.Script import Script

# ----------------------------------------------------------------------
class TestApp(Script):
  """
  Test application.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="testapp"):
    """
    Constructor.
    """
    Script.__init__(self, name)
    return


  def main(self, *args, **kwds):
    """
    Run the application.
    """
    from pylith.utils.PetscManager import PetscManager
    manager = PetscManager()
    manager.initialize()

    self._runTest()

    manager.finalize()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _runTest(self):
    """
    Run test.
    """

    #filenameIn = "data/tet4.exo"
    filenameIn = "tri3_200m_gradient.exo"

    from pylith.perf.MemoryLogger import MemoryLogger
    self.logger = MemoryLogger()
    self.logger._configure()

    from pylith.utils.petsc import MemoryLogger
    sieveLogger =  MemoryLogger.singleton()

    from pylith.topology.topology import MeshOps_numMaterialCells

    from pylith.mpi.Communicator import petsc_comm_world
    comm = petsc_comm_world()

    self._showStatus("Before reading mesh")

    # ------------------------------------------------------------
    # Read mesh
    from pylith.meshio.MeshIOCubit import MeshIOCubit
    io = MeshIOCubit()
    io.inventory.filename = filenameIn
    io.inventory.useNames = False
    io._configure()

    from pylith.materials.ElasticIsotropic3D import ElasticIsotropic3D
    material = ElasticIsotropic3D()
    material.inventory.label = "Elastic material"
    material.inventory.id = 1
    material._configure()

    mesh = io.read(debug=False, interpolate=False)

    self.logger.logMesh(mesh.memLoggingStage, mesh)
    material.ncells = MeshOps_numMaterialCells(mesh, material.id())
    self.logger.logMaterial(mesh.memLoggingStage, material)
    

    self._showStatus("After reading mesh")    

    # ------------------------------------------------------------
    # Reorder mesh
    from pylith.topology.ReverseCuthillMcKee import ReverseCuthillMcKee
    ordering = ReverseCuthillMcKee()
    ordering.reorder(mesh)

    # Expect no memory allocation
    self.logger.memory["MeshReordering"] = 0

    self._showStatus("After reordering mesh")


    # ------------------------------------------------------------
    # Distribute mesh
    if comm.size > 1:
      from pylith.topology.Distributor import Distributor
      from spatialdata.units.Nondimensional import Nondimensional
      distributor = Distributor()
      normalizer = Nondimensional()
      dmesh = distributor.distribute(mesh, normalizer)
      dmesh.memLoggingStage = "DistributedMesh"

      self.logger.logMesh(mesh.memLoggingStage, mesh)
      material.ncells = MeshOps_numMaterialCells(mesh, material.id())
      self.logger.logMaterial(mesh.memLoggingStage, material)
    
      self.logger.logMesh(dmesh.memLoggingStage, dmesh)
      material.ncells = MeshOps_numMaterialCells(dmesh, material.id())
      self.logger.logMaterial(dmesh.memLoggingStage, material)

      self._showStatus("After distributing mesh")

      mesh = dmesh
      
      
    # Refine mesh (if necessary)
    from pylith.topology.RefineUniform import RefineUniform
    refiner = RefineUniform()
    rmesh = refiner.refine(mesh)
    rmesh.memLoggingStage = "RefinedMesh"

    print "Unrefined mesh logging stage",mesh.memLoggingStage
    self.logger.logMesh(mesh.memLoggingStage, mesh)
    material.ncells = MeshOps_numMaterialCells(mesh, material.id())
    self.logger.logMaterial(mesh.memLoggingStage, material)
    
    self.logger.logMesh(rmesh.memLoggingStage, rmesh)
    material.ncells = MeshOps_numMaterialCells(rmesh, material.id())
    self.logger.logMaterial(rmesh.memLoggingStage, material)

    self._showStatus("After refining mesh")


    return


  def _showStatus(self, stage):
    from pylith.utils.profiling import resourceUsageString
    from pylith.mpi.Communicator import petsc_comm_world
    comm = petsc_comm_world()

    if comm.rank == 0:
      print "\n----------------------------------------------------------------------"
      print "STAGE: %s" % stage
      print "----------------------------------------------------------------------"
    for irank in xrange(comm.size):
      comm.barrier()
      if comm.rank == irank:
        print "\nPROCESSOR %d" % comm.rank
        print "\nStatus from ps: %s\n" % resourceUsageString()
        self.logger.show()

      comm.barrier()
    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
