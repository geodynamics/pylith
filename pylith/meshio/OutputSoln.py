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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pyre/meshio/OutputSoln.py
##
## @brief Python object for managing output of finite-element
## solution information.
##
## Factory: output_manager

from OutputManager import OutputManager

# OutputSoln class
class OutputSoln(OutputManager):
  """
  Python object for managing output of finite-element solution
  information.

  @class Inventory
  Python object for managing OutputSoln facilities and properties.
  
  \b Properties
  @li \b vertex_data_fields Names of vertex data fields to output.
  @li \b cell_info_fields Names of cell info fields to output.
  
  \b Facilities
  @li None

  Factory: mesh_output_manager
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  vertexDataFields = pyre.inventory.list("vertex_data_fields", 
                                         default=["displacement"])
  vertexDataFields.meta['tip'] = "Names of vertex data fields to output."
  
  cellInfoFields = pyre.inventory.list("cell_info_fields", default=[])
  cellInfoFields.meta['tip'] = "Names of cell info fields to output."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="outputsoln"):
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

    import weakref
    self.mesh = weakref.ref(mesh)
    OutputManager.initialize(self, normalizer)

    self._eventLogger.eventEnd(logEvent)
    return


  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.mesh(), None, None)


  def getVertexField(self, name, fields):
    """
    Get vertex field.
    """
    # :TODO: Clean this up for multiple fields

    buffer = None
    if name == "displacement":
      field = fields.get("disp(t)")
      if not fields.hasField("buffer (vector)"):
        fields.add("buffer (vector)", "buffer")
      buffer = fields.get("buffer (vector)")
      buffer.copySubfield(field, "displacement")
    elif name == "velocity":
      field = fields.get("velocity(t)")
      if not fields.hasField("buffer (vector)"):
        fields.add("buffer (vector)", "buffer")
      buffer = fields.get("buffer (vector)")
      buffer.copySubfield(field, "displacement")
      buffer.label(field.label()) # :KLUDGE: Fix for multiple fields
      buffer.scale(field.scale()) # :KLUDGE: Fix for multiple fields
    else:
      raise ValueError, "Vertex field '%s' not available." % name

    buffer.dimensionalizeOkay(True)
    return buffer


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    OutputManager._configure(self)
    self.vertexDataFields = self.inventory.vertexDataFields
    self.cellInfoFields = self.inventory.cellInfoFields
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputSoln()


# End of file 
