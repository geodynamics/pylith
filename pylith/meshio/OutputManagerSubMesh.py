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

## @file pyre/meshio/OutputManagerSubMesh.py
##
## @brief Python abstract base class for managing output of
## finite-element information.
##
## Factory: output_manager

from OutputManager import OutputManager
from meshio import SubMeshOutputManager as ModuleOutputManager

# OutputManagerSubMsh class
class OutputManagerSubMesh(OutputManager, ModuleOutputManager):
  """
  Python abstract base class for managing output of finite-element
  information.

  \b Properties
  @li None
  
  \b Facilities
  @li \b writer Writer for data.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  from DataWriterVTKSubSubMesh import DataWriterVTKSubSubMesh
  writer = pyre.inventory.facility("writer", factory=DataWriterVTKSubSubMesh,
                                   family="data_writer")
  writer.meta['tip'] = "Writer for data."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputmanagersubmesh"):
    """
    Constructor.
    """
    OutputManager.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    OutputManager._configure(self)
    ModuleOutputManager.coordsys(self, self.inventory.coordsys)
    ModuleOutputManager.writer(self, self.inventory.writer)
    if None != self.vertexFilter.filter:
      ModuleOutputManager.vertexFilter(self, self.inventory.vertexFilter)
    if None != self.cellFilter.filter:
      ModuleOutputManager.cellFilter(self, self.inventory.cellFilter)
    return


  def _createModuleObj(self):
    """
    Create handle to C++ object.
    """
    ModuleOutputManager.__init__(self)
    return


  def _open(self, mesh, nsteps, label, labelId):
    """
    Call C++ open();
    """
    if label != None and labelId != None:
      ModuleOutputManager.open(self, mesh, nsteps, label, labelId)
    else:
      ModuleOutputManager.open(self, mesh, nsteps)
    return


  def _openTimeStep(self, t, mesh, label, labelId):
    """
    Call C++ openTimeStep();
    """
    if label != None and labelId != None:
      ModuleOutputManager.openTimeStep(self, t, mesh, label, labelId)
    else:
      ModuleOutputManager.openTimeStep(self, t, mesh)
    return


  def _appendVertexField(self, t, field):
    """
    Call C++ appendVertexField();
    """
    ModuleOutputManager.appendVertexField(self, t, field)
    return

  def _appendCellField(self, t, field):
    """
    Call C++ appendCellField();
    """
    ModuleOutputManager.appendCellField(self, t, field)
    return


  def _closeTimeStep(self):
    """
    Call C++ closeTimeStep().
    """
    ModuleOutputManager.closeTimeStep(self)
    return


  def _close(self):
    """
    Call C++ close().
    """
    ModuleOutputManager.close(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputManagerSubMesh()


# End of file 
