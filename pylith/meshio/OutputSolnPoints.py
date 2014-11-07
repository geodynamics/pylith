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
# Copyright (c) 2010-2014 University of California, Davis
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

  Factory: output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  vertexDataFields = pyre.inventory.list("vertex_data_fields", default=["displacement"])
  vertexDataFields.meta['tip'] = "Names of vertex data fields to output."
  
  from PointsList import PointsList
  reader = pyre.inventory.facility("reader", factory=PointsList, family="points_list")
  reader.meta['tip'] = "Reader for points list."


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
  

  def initialize(self, mesh, normalizer):
    """
    Initialize output manager.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)    

    OutputManager.initialize(self, normalizer)

    # Read points
    self.stations,points = self.reader.read()
    
    # Convert to mesh coordinate system
    from spatialdata.geocoords.Converter import convert
    convert(points, mesh.coordsys(), self.coordsys)

    ModuleOutputSolnPoints.setupInterpolator(self, mesh, points, normalizer)
    self.mesh = ModuleOutputSolnPoints.pointsMesh(self)

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

    ModuleOutputSolnPoints.writePointNames(self, self.stations)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputSolnPoints()


# End of file 
