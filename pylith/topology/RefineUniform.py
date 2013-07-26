#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/topology/RefineUniform.py
##
## @brief Python manager for uniform global refinement of mesh in
## parallel.
##
## Factory: mesh_refiner.

from MeshRefiner import MeshRefiner
from topology import RefineUniform as ModuleRefineUniform

# RefineUniform class
class RefineUniform(MeshRefiner, ModuleRefineUniform):
  """
  Python manager for uniform global refinement of mesh in parallel.

  Factory: mesh_refiner
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  levels = pyre.inventory.int("levels", default=1, validator=pyre.inventory.greaterEqual(1))
  levels.meta['tip'] = "Number of refinement levels."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="refineuniform"):
    """
    Constructor.
    """
    MeshRefiner.__init__(self, name)
    self._createModuleObj()
    return


  def refine(self, mesh):
    """
    Refine mesh.
    """
    self._setupLogging()
    logEvent = "%srefine" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    newMesh = mesh
    from Mesh import Mesh
    newMesh = Mesh(dim=mesh.dimension(), comm=mesh.comm())
    newMesh.debug(mesh.debug())
    newMesh.coordsys(mesh.coordsys())
    ModuleRefineUniform.refine(self, newMesh, mesh, self.levels)
    if not newMesh == mesh:
      #from pylith.utils.petsc import MemoryLogger
      #memoryLogger =  MemoryLogger.singleton()
    
      #memoryLogger.stagePush(mesh.memLoggingStage)      
      mesh.cleanup()
      #memoryLogger.stagePop()

    self._eventLogger.eventEnd(logEvent)
    return newMesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    MeshRefiner._configure(self)
    self.levels = self.inventory.levels
    return


  def _createModuleObj(self):
    """
    Create handle to C++ object.
    """
    ModuleRefineUniform.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_refiner():
  """
  Factory associated with RefineUniform.
  """
  return RefineUniform()


# End of file 
