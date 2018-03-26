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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/faults/FaultCohesive.py
##

## @brief Python abstract base class for a fault surface implemented
## with cohesive elements.
##
## Factory: fault

from pylith.utils.PetscComponent import PetscComponent
from faults import FaultCohesive as ModuleFaultCohesive

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


# FaultCohesive class
class FaultCohesive(Fault, ModuleFaultCohesive):
  """
  Python abstract base class for a fault surface implemeted with
  cohesive elements.

  Inventory

  @class Inventory
  Python object for managing FaultCohesive facilities and properties.
  
  \b Properties
  @li \b id Fault identifier
  @li \b label Label identifier for fault.
  @li \b edge Label identifier for buried fault edges.
  @li \b up_dir Up-dip or up direction
    (perpendicular to along-strike and not collinear with fault normal;
    applies to fault surfaces in 2-D and 3-D).
  
  \b Facilities

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
  faultEdge.meta['tip'] = "Label identifier for buried fault edges."
  
  upDir = pyre.inventory.list("up_dir", default=[0.0, 0.0, 1.0], validator=validateDir)
  upDir.meta['tip'] = "Up-dip or up direction " \
      "(perpendicular to along-strike and not collinear " \
      "with fault normal; applies to fault surfaces " \
      "in 2-D and 3-D)."


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


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    PetscComponent._configure(self)
    self.faultQuadrature = self.inventory.faultQuadrature
    self.upDir = map(float, self.inventory.upDir)
    ModuleFault.id(self, self.inventory.matId)
    ModuleFault.label(self, self.inventory.faultLabel)
    ModuleFault.edge(self, self.inventory.faultEdge)
    self.perfLogger = self.inventory.perfLogger
    return

  
  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    raise NotImplementedError("Please implement _createModuleObj() in derived class.")
  
  
# End of file 
