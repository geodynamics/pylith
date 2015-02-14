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

## @file pylith/faults/Fault.py
##

## @brief Python abstract base class for a fault surface.
##
## This implementation of a fault associates both physical
## properties and a quadrature scheme with the fault.
##
## Factory: fault

from pylith.utils.PetscComponent import PetscComponent
from faults import Fault as ModuleFault

# Validator for label
def validateLabel(value):
  """
  Validate label for group/nodeset/pset.
  """
  if 0 == len(value):
    raise ValueError("Label for group/nodeset/pset in mesh not specified.")
  return value


# Validator for direction
def validateDir(value):
  """
  Validate direction.
  """
  msg = "Direction must be a 3 component vector (list)."
  if not isinstance(value, list):
    raise ValueError(msg)
  if 3 != len(value):
    raise ValueError(msg)
  try:
    nums = map(float, value)
  except:
    raise ValueError(msg)
  return nums


# Fault class
class Fault(PetscComponent, ModuleFault):
  """
  Python abstract base class for a fault surface.

  This implementation of a fault associates both physical
  properties and a quadrature scheme with the fault.

  Inventory

  \b Properties
  @li \b id Fault identifier
  @li \b label Label identifier for fault.
  @li \b up_dir Up-dip or up direction
    (perpendicular to along-strike and not collinear with fault normal;
    applies to fault surfaces in 2-D and 3-D).
  
  \b Facilities
  @li \b quadrature Quadrature object for numerical integration

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory
  
  matId = pyre.inventory.int("id", default=100)
  matId.meta['tip'] = "Fault identifier (must be unique across all faults " \
      "and materials)."
  
  faultLabel = pyre.inventory.str("label", default="", validator=validateLabel)
  faultLabel.meta['tip'] = "Label identifier for fault."
  
  faultEdge = pyre.inventory.str("edge", default="")
  faultEdge.meta['tip'] = "Label identifier for fault edge."
  
  upDir = pyre.inventory.list("up_dir", default=[0, 0, 1],
                              validator=validateDir)
  upDir.meta['tip'] = "Up-dip or up direction " \
      "(perpendicular to along-strike and not collinear " \
      "with fault normal; applies to fault surfaces " \
      "in 2-D and 3-D)."
  
  from pylith.feassemble.Quadrature import Quadrature
  faultQuadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
  faultQuadrature.meta['tip'] = "Quadrature object for numerical integration."
  
  from pylith.perf.MemoryLogger import MemoryLogger
  perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger",
                                       factory=MemoryLogger)
  perfLogger.meta['tip'] = "Performance and memory logging."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fault"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="fault")
    self._createModuleObj()
    self.mesh = None
    self.output = None
    return


  def preinitialize(self, mesh):
    """
    Setup fault.
    """
    import weakref
    self.mesh = weakref.ref(mesh)
    
    self.faultQuadrature.preinitialize(mesh.coordsys().spaceDim())

    if None != self.output:
      self.output.preinitialize(self)

    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    faultDim = self.mesh().dimension() - 1
    if faultDim != self.faultQuadrature.cell.cellDim:
      raise ValueError, \
            "Quadrature is incompatible with fault surface.\n" \
            "Dimensions for quadrature: %d, dimensions of fault: %d" % \
            (self.faultQuadrature.cell.cellDim, faultDim)

    if None != self.output:
      self.output.verifyConfiguration(self.mesh())

    self._eventLogger.eventEnd(logEvent)
    return
  

  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize fault.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    self.faultQuadrature.initialize()
    ModuleFault.initialize(self, self.mesh(), self.upDir)

    if None != self.output:
      self.output.initialize(normalizer, self.faultQuadrature)
      self.output.writeInfo()
      self.output.open(totalTime, numTimeSteps)

    self._eventLogger.eventEnd(logEvent)
    return


  def writeData(self, t, fields):
    """
    Write data at time t.
    """
    logEvent = "%swrite" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()

    if 0 == comm.rank:
      self._info.log("Writing fault data.")
    self.output.writeData(t, fields)

    self._eventLogger.eventEnd(logEvent)
    return


  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.faultMesh(), None, None)


  def getVertexField(self, name, fields=None):
    """
    Get vertex field.
    """
    raise NotImplementedError("Fault.getVertexField() not implemented.")
    return


  def getCellField(self, name, fields=None):
    """
    Get cell field.
    """
    raise NotImplementedError("Fault.getCellField() not implemented.")
    return


  def finalize(self):
    """
    Cleanup.
    """
    self._modelMemoryUse()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    try:
      PetscComponent._configure(self)
      self.faultQuadrature = self.inventory.faultQuadrature
      self.upDir = map(float, self.inventory.upDir)
      ModuleFault.id(self, self.inventory.matId)
      ModuleFault.label(self, self.inventory.faultLabel)
      ModuleFault.edge(self, self.inventory.faultEdge)
      self.perfLogger = self.inventory.perfLogger
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring fault "
                       "(%s):\n%s" % (aliases, err.message))
    return

  
  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    raise NotImplementedError("Please implement _createModuleObj() in " \
                              "derived class.")
  
  
  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    raise NotImplementedError, \
          "Please implement _modelModelUse() in derived class."
    return


# End of file 
