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

## @file pyre/meshio/OutputSolnSubset.py
##
## @brief Python object for managing output of finite-element solution
## information over a subdomain.
##
## Factory: output_manager

from OutputManager import OutputManager
from meshio import OutputSolnSubset as ModuleOutputSolnSubset

# Validator for label
def validateLabel(value):
  """
  Validate label for group/nodeset/pset.
  """
  if 0 == len(value):
    raise ValueError("Label for group/nodeset/pset in mesh not specified.")
  return value


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
  
  label = pyre.inventory.str("label", default="", validator=validateLabel)
  label.meta['tip'] = "Label identifier for subdomain."

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
    try:
      OutputManager._configure(self)
      ModuleOutputSolnSubset.label(self, self.label)
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring output over boundary "
                       "(%s):\n%s" % (aliases, err.message))

    return


  def _createModuleObj(self):
    """
    Create handle to C++ object.
    """
    ModuleOutputSolnSubset.__init__(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
  """
  Factory associated with OutputManager.
  """
  return OutputSolnSubset()


# End of file 
