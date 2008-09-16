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

# RefineUniform class
class RefineUniform(MeshRefiner):
  """
  Python manager for uniform global refinement of mesh in parallel.

  Factory: mesh_refiner
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MeshRefiner.Inventory):
    """
    Python object for managing RefineUniform facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing RefineUniform facilities and properties.
    ##
    ## \b Properties
    ## @li \b levels Number of refinement levels.
    ##
    ## \b Facilities
    ## @li None

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
    return


  def refine(self, mesh):
    """
    Refine mesh.
    """
    self._setupLogging()
    logEvent = "%srefine" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._createCppHandle()
    
    from Mesh import Mesh
    newMesh = Mesh()
    assert(None != self.cppHandle)
    newMesh.cppHandle = self.cppHandle.refine(mesh.cppHandle,
                                              self.levels)
    newMesh.coordsys = mesh.coordsys

    if self.debug:
      self.dataWriter.initialize()
      self.cppHandle.write(self.dataWriter.cppHandle,
                           newMesh.cppHandle, newMesh.coordsys.cppHandle)

    self._logger.eventEnd(logEvent)
    return newMesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.levels = self.inventory.levels
    return


  def _createCppHandle(self):
    """
    Create handle to C++ object.
    """
    if None == self.cppHandle:
      import pylith.topology.topology as bindings
      self.cppHandle = bindings.RefineUniform()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_distributor():
  """
  Factory associated with RefineUniform.
  """
  return RefineUniform()


# End of file 
