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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
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
    from pylith.utils.NullComponent import NullComponent

    OutputManager._configure(self)
    ModuleOutputManager.coordsys(self, self.inventory.coordsys)
    ModuleOutputManager.writer(self, self.inventory.writer)
    if not isinstance(self.inventory.vertexFilter, NullComponent):
      ModuleOutputManager.vertexFilter(self, self.inventory.vertexFilter)
    if not isinstance(self.inventory.cellFilter, NullComponent):
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


  def _appendVertexField(self, t, field, mesh):
    """
    Call C++ appendVertexField();
    """
    ModuleOutputManager.appendVertexField(self, t, field, mesh)
    return


  def _appendCellField(self, t, field, label, labelId):
    """
    Call C++ appendCellField();
    """
    if label != None and labelId != None:
      ModuleOutputManager.appendCellField(self, t, field, label, labelId)
    else:
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
