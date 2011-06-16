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

## @file pyre/meshio/OutputSolnSubset.py
##
## @brief Python object for managing output of finite-element solution
## information over a subdomain.
##
## Factory: output_manager

from OutputManager import OutputManager
from meshio import OutputSolnSubset as ModuleOutputSolnSubset

# OutputSolnSubset class
class OutputSolnSubset(OutputManager, ModuleOutputSolnSubset):
  """
  Python object for managing output of finite-element solution
  information over a subdomain.

  @class Inventory
  Python object for managing OutputSolnSubset facilities and properties.
  
  \b Properties
  @li \b vertex_data_fields Names of vertex data fields to output.
  @li \b label Name identifier for subdomain.
  
  \b Facilities
  @li \b writer Writer for data.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  vertexDataFields = pyre.inventory.list("vertex_data_fields", 
                                         default=["displacement"])
  vertexDataFields.meta['tip'] = "Names of vertex data fields to output."
  
  label = pyre.inventory.str("label", default="")
  label.meta['tip'] = "Label identifier for subdomain."

  from DataWriterVTKSubMesh import DataWriterVTKSubMesh
  writer = pyre.inventory.facility("writer", factory=DataWriterVTKSubMesh,
                                 family="data_writer")
  writer.meta['tip'] = "Writer for data."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputsolnsubset"):
    """
    Constructor.
    """
    OutputManager.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': ["displacement","velocity"]},
         'cell': \
           {'info': [],
            'data': []}}
    return


  def preinitialize(self):
    """
    Do
    """
    OutputManager.preinitialize(self, dataProvider=self)
    return
  

  def verifyConfiguration(self, mesh):
    """
    Verify compatibility of configuration.
    """
    OutputManager.verifyConfiguration(self, mesh)
    ModuleOutputSolnSubset.verifyConfiguration(self, mesh)
    return


  def initialize(self, mesh, normalizer):
    """
    Initialize output manager.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)    

    self.submesh = self.subdomainMesh(mesh)
    OutputManager.initialize(self, normalizer)

    self._eventLogger.eventEnd(logEvent)
    return


  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.submesh, None, None)


  def getVertexField(self, name, fields):
    """
    Get vertex field.
    """
    field = None
    fieldType = None
    if name == "displacement":
      field = fields.get("disp(t)")
    elif name == "velocity":
      field = fields.get("velocity(t)")
    else:
      raise ValueError, "Vertex field '%s' not available." % name
    return field


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    OutputManager._configure(self)
    ModuleOutputSolnSubset.label(self, self.label)
    ModuleOutputSolnSubset.coordsys(self, self.inventory.coordsys)
    ModuleOutputSolnSubset.writer(self, self.inventory.writer)
    from pylith.utils.NullComponent import NullComponent
    if not isinstance(self.inventory.vertexFilter, NullComponent):
      ModuleOutputSolnSubset.vertexFilter(self, self.inventory.vertexFilter)
    if not isinstance(self.inventory.cellFilter, NullComponent):
      ModuleOutputSolnSubset.cellFilter(self, self.inventory.cellFilter)
    return


  def _createModuleObj(self):
    """
    Create handle to C++ object.
    """
    ModuleOutputSolnSubset.__init__(self)
    return


  def _open(self, mesh, nsteps, label, labelId):
    """
    Call C++ open();
    """
    if label != None and labelId != None:
      ModuleOutputSolnSubset.open(self, mesh, nsteps, label, labelId)
    else:
      ModuleOutputSolnSubset.open(self, mesh, nsteps)
    return


  def _openTimeStep(self, t, mesh, label, labelId):
    """
    Call C++ openTimeStep();
    """
    if label != None and labelId != None:
      ModuleOutputSolnSubset.openTimeStep(self, t, mesh, label, labelId)
    else:
      ModuleOutputSolnSubset.openTimeStep(self, t, mesh)
    return


  def _appendVertexField(self, t, field, mesh):
    """
    Call C++ appendVertexField();
    """
    ModuleOutputSolnSubset.appendVertexField(self, t, field, mesh)
    return

  def _appendCellField(self, t, field):
    """
    Call C++ appendCellField();
    """
    ModuleOutputSolnSubset.appendCellField(self, t, field)
    return


  def _closeTimeStep(self):
    """
    Call C++ closeTimeStep().
    """
    ModuleOutputSolnSubset.closeTimeStep(self)
    return


  def _close(self):
    """
    Call C++ close().
    """
    ModuleOutputSolnSubset.close(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputSolnSubset()


# End of file 
