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

## @file pyre/meshio/OutputSolnPoints.py
##
## @brief Python object for managing output of finite-element solution
## information over a subdomain.
##
## Factory: output_manager

from OutputManager import OutputManager
from meshio import OutputSolnPoints as ModuleOutputSolnPoints

# Validator for filename
def validateFilename(value):
  """
  Validate filename with list of points.
  """
  if 0 == len(value):
    raise ValueError("Filename for list of points not specified.")
  return value


# OutputSolnPoints class
class OutputSolnPoints(OutputManager, ModuleOutputSolnPoints):
  """
  Python object for managing output of finite-element solution
  information over a subdomain.

  @class Inventory
  Python object for managing OutputSolnPoints facilities and properties.
  
  \b Properties
  @li \b vertex_data_fields Names of vertex data fields to output.
  
  \b Facilities
  @li \b reader Reader for list of points.
  @li \b writer Writer for data.

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  vertexDataFields = pyre.inventory.list("vertex_data_fields", 
                                         default=["displacement"])
  vertexDataFields.meta['tip'] = "Names of vertex data fields to output."
  
  from PointsList import PointsList
  reader = pyre.inventory.facility("reader", factory=PointsList, family="points_list")
  reader.meta['tip'] = "Reader for points list."

  from DataWriterVTKPoints import DataWriterVTKPoints
  writer = pyre.inventory.facility("writer", factory=DataWriterVTKPoints,
                                 family="data_writer")
  writer.meta['tip'] = "Writer for data."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputsolnpoints"):
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
    ModuleOutputSolnPoints.verifyConfiguration(self, mesh)
    return


  def initialize(self, mesh, normalizer):
    """
    Initialize output manager.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)    

    points = self.reader.read()
    ModuleOutputSolnPoints.setupInterpolator(self, mesh, points)
    self.mesh = ModuleOutputSolnPoints.createPointsMesh(self)
    OutputManager.initialize(self, normalizer)

    self._eventLogger.eventEnd(logEvent)
    return


  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.mesh, None, None)


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
    try:
      OutputManager._configure(self)
      ModuleOutputSolnPoints.label(self, self.label)
      ModuleOutputSolnPoints.coordsys(self, self.inventory.coordsys)
      ModuleOutputSolnPoints.writer(self, self.inventory.writer)
      from pylith.utils.NullComponent import NullComponent
      if not isinstance(self.inventory.vertexFilter, NullComponent):
        ModuleOutputSolnPoints.vertexFilter(self, self.inventory.vertexFilter)
      if not isinstance(self.inventory.cellFilter, NullComponent):
        ModuleOutputSolnPoints.cellFilter(self, self.inventory.cellFilter)
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring output over points "
                       "(%s):\n%s" % (aliases, err.message))

    return


  def _createModuleObj(self):
    """
    Create handle to C++ object.
    """
    ModuleOutputSolnPoints.__init__(self)
    return


  def _open(self, mesh, nsteps, label, labelId):
    """
    Call C++ open();
    """
    if label != None and labelId != None:
      ModuleOutputSolnPoints.open(self, mesh, nsteps, label, labelId)
    else:
      ModuleOutputSolnPoints.open(self, mesh, nsteps)
    return


  def _openTimeStep(self, t, mesh, label, labelId):
    """
    Call C++ openTimeStep();
    """
    if label != None and labelId != None:
      ModuleOutputSolnPoints.openTimeStep(self, t, mesh, label, labelId)
    else:
      ModuleOutputSolnPoints.openTimeStep(self, t, mesh)
    return


  def _appendVertexField(self, t, field, mesh):
    """
    Call C++ appendVertexField();
    """
    ModuleOutputSolnPoints.appendVertexField(self, t, field, mesh)
    return

  def _appendCellField(self, t, field):
    """
    Call C++ appendCellField();
    """
    raise NotImplementedError("Output of cell field not implemented for arbitrary points.")
    return


  def _closeTimeStep(self):
    """
    Call C++ closeTimeStep().
    """
    ModuleOutputSolnPoints.closeTimeStep(self)
    return


  def _close(self):
    """
    Call C++ close().
    """
    ModuleOutputSolnPoints.close(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputSolnPoints()


# End of file 
