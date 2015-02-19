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

## @file pylith/bc/DirichletBoundary.py
##
## @brief Python object for managing a Dirichlet (prescribed
## displacements) boundary condition with points on a surface.
##
## Factory: boundary_condition

from DirichletBC import DirichletBC
from bc import DirichletBoundary as ModuleDirichletBoundary

from pylith.utils.NullComponent import NullComponent

# DirichletBoundary class
class DirichletBoundary(DirichletBC, ModuleDirichletBoundary):
  """
  Python object for managing a DirichletBoundary (prescribed displacements)
  boundary condition with points on a surface.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(DirichletBC.Inventory):
    """
    Python object for managing DirichletBoundary facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing DirichletBoundary facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b output Output manager associated with diagnostic output.

    import pyre.inventory

    from pylith.meshio.OutputDirichlet import OutputDirichlet
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputDirichlet)
    output.meta['tip'] = "Output manager associated with diagnostic output."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="dirichletboundary"):
    """
    Constructor.
    """
    DirichletBC.__init__(self, name)
    self._loggingPrefix = "DiBC "
    self.availableFields = \
        {'vertex': \
           {'info': [],
            'data': []},
         'cell': \
           {'info': [],
            'data': []}}
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    DirichletBC.preinitialize(self, mesh)
    self.output.preinitialize(self)

    fields = []
    if not isinstance(self.inventory.dbInitial, NullComponent):
      fields += ["initial_value"]
    if not isinstance(self.inventory.dbRate, NullComponent):
      fields += ["rate_of_change", "rate_start_time"]
    if not isinstance(self.inventory.dbChange, NullComponent):
      fields += ["change_in_value", "change_start_time"]
    self.availableFields['vertex']['info'] += fields
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    DirichletBC.verifyConfiguration(self)
    self.output.verifyConfiguration(self.mesh())

    self._eventLogger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize DirichletBoundary boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Initializing Dirichlet boundary '%s'." % self.label())

    DirichletBC.initialize(self, totalTime, numTimeSteps, normalizer)
    self.output.initialize(normalizer)
    self.output.writeInfo()

    self._eventLogger.eventEnd(logEvent)    
    return
  

  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    label = ""
    labelId = 0
    return (self.boundaryMesh(), label, labelId)


  def getVertexField(self, name, fields=None):
    """
    Get vertex field.
    """
    if None == fields:
      field = self.vertexField(name)
    else:
      field = self.vertexField(name, fields)
    return field


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    DirichletBC._configure(self)
    self.output = self.inventory.output
    return


  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    ModuleDirichletBoundary.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with DirichletBoundary.
  """
  return DirichletBoundary()

  
# End of file 
