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

## @file pylith/bc/Neumann.py
##
## @brief Python object for managing traction boundary conditions.
##
## Factory: boundary_condition

from BoundaryCondition import BoundaryCondition
from TimeDependent import TimeDependent
from pylith.feassemble.Integrator import Integrator
from bc import Neumann as ModuleNeumann

from pylith.utils.NullComponent import NullComponent

# Neumann class
class Neumann(BoundaryCondition, 
              TimeDependent,
              Integrator, 
              ModuleNeumann):
  """
  Python object for managing traction boundary conditions.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory
  
  from pylith.feassemble.Quadrature import Quadrature
  bcQuadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
  bcQuadrature.meta['tip'] = "Quadrature object for numerical integration."
  
  from pylith.meshio.OutputNeumann import OutputNeumann
  output = pyre.inventory.facility("output", family="output_manager",
                                   factory=OutputNeumann)
  output.meta['tip'] = "Output manager associated with diagnostic output."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="neumann"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    Integrator.__init__(self)
    TimeDependent.__init__(self)
    self._loggingPrefix = "NeBC "
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
    BoundaryCondition.preinitialize(self, mesh)
    Integrator.preinitialize(self, mesh)
    self.bcQuadrature.preinitialize(mesh.coordsys().spaceDim())
    self.quadrature(self.bcQuadrature)
    self.createSubMesh(mesh)
    self.output.preinitialize(self)

    fields = []
    if not isinstance(self.inventory.dbInitial, NullComponent):
      fields += ["initial_value"]
    if not isinstance(self.inventory.dbRate, NullComponent):
      fields += ["rate_of_change", "rate_start_time"]
    if not isinstance(self.inventory.dbChange, NullComponent):
      fields += ["change_in_value", "change_start_time"]
    self.availableFields['cell']['info'] += fields
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    BoundaryCondition.verifyConfiguration(self, self.mesh())
    Integrator.verifyConfiguration(self)
    if self.bcQuadrature.cellDim() != self.mesh().dimension()-1:
        raise ValueError, \
              "Quadrature scheme and mesh are incompatible.\n" \
              "Dimension for quadrature: %d\n" \
              "Dimension of mesh boundary '%s': %d" % \
              (self.bcQuadrature.cellDim(),
               self.label(), self.mesh().dimension()-1)    
    self.output.verifyConfiguration(self.mesh())
    ModuleNeumann.verifyConfiguration(self, self.mesh())

    self._eventLogger.eventEnd(logEvent)
    return
  

  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize Neumann boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)
    
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Initializing Neumann boundary '%s'." % self.label())

    Integrator.initialize(self, totalTime, numTimeSteps, normalizer)
    BoundaryCondition.initialize(self, totalTime, numTimeSteps, normalizer)

    self.output.initialize(normalizer, self.bcQuadrature)
    self.output.writeInfo()

    self._eventLogger.eventEnd(logEvent)
    return
  

  def finalize(self):
    """
    Cleanup.
    """
    BoundaryCondition.finalize(self)
    Integrator.finalize(self)
    self._modelMemoryUse()
    return
  

  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.boundaryMesh(), None, None)


  def getCellField(self, name, fields=None):
    """
    Get vertex field.
    """
    if None == fields:
      field = self.cellField(name)
    else:
      field = self.cellField(name, fields)
    return field


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    TimeDependent._configure(self)
    self.bcQuadrature = self.inventory.bcQuadrature
    self.output = self.inventory.output
    return


  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    ModuleNeumann.__init__(self)
    return
  

  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    self.perfLogger.logFields("BoundaryConditions", self.parameterFields())
    return


# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with Neumann.
  """
  return Neumann()

  
# End of file 
