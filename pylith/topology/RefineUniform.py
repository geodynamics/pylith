#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
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

  levels = pyre.inventory.int("levels", default=2,
                              validator=pyre.inventory.choice([2, 4]))
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

    from Mesh import Mesh
    newMesh = Mesh(mesh.dimension(), mesh.getComm())
    newMesh.debug(mesh.debug())
    newMesh.coordsys(mesh.coordsys())
    ModuleRefineUniform.refine(self, newMesh, mesh, self.levels)

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
