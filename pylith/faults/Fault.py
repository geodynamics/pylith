#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
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
  @li \b name Name of fault
  @li \b up_dir Up-dip or up direction
    (perpendicular to along-strike and not collinear with fault normal;
    only applies to fault surfaces in a 3-D domain).
  @li \b normal_dir General preferred direction for fault normal
    (used to pick which of two possible normal directions for
    interface; only applies to fault surfaces in a 3-D domain).
  @li \b mat_db Spatial database for bulk material properties
    (used in improving conditioning of Jacobian matrix).
  
  \b Facilities
  @li \b quadrature Quadrature object for numerical integration

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory
  
  id = pyre.inventory.int("id", default=100)
  id.meta['tip'] = "Fault identifier (must be unique across all faults " \
      "and materials)."
  
  label = pyre.inventory.str("label", default="")
  label.meta['tip'] = "Name of fault."
  
  upDir = pyre.inventory.list("up_dir", default=[0, 0, 1],
                              validator=validateDir)
  upDir.meta['tip'] = "Up-dip or up direction " \
      "(perpendicular to along-strike and not collinear " \
      "with fault normal; only applies to fault surface " \
      "in a 3-D domain)."
  
  normalDir = pyre.inventory.list("normal_dir", default=[1, 0, 0],
                                  validator=validateDir)
  normalDir.meta['tip'] = "General preferred direction for fault normal " \
      "(used to pick which of two possible normal directions for " \
      "interface; only applies to fault surfaces in a 3-D domain)."
  
  from pylith.feassemble.Quadrature import SubMeshQuadrature
  quadrature = pyre.inventory.facility("quadrature", factory=SubMeshQuadrature)
  quadrature.meta['tip'] = "Quadrature object for numerical integration."
  
  from spatialdata.spatialdb.SimpleDB import SimpleDB
  matDB = pyre.inventory.facility("mat_db", family="spatial_database",
                                  factory=SimpleDB)
  matDB.meta['tip'] = "Spatial database for bulk material properties " \
      "(used in improving conditioning of Jacobian matrix)."

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
    self.mesh = mesh
    
    self.faultQuadrature.preinitialize(mesh.coordsys().spaceDim())

    if None != self.output:
      self.output.preinitialize(self)

    return
  

  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    faultDim = self.mesh.dimension() - 1
    if faultDim != self.faultQuadrature.cell.cellDim:
      raise ValueError, \
            "Quadrature is incompatible with fault surface.\n" \
            "Dimensions for quadrature: %d, dimensions of fault: %d" % \
            (self.faultQuadrature.cell.cellDim, faultDim)

    if None != self.output:
      self.output.verifyConfiguration(self.mesh)

    self._logger.eventEnd(logEvent)
    return
  

  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize fault.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self.faultQuadrature.initialize()
    ModuleFault.initialize(self, 
                           self.mesh, self.upDir, self.normalDir, self.matDB)

    if None != self.output:
      self.output.initialize(normalizer, self.faultQuadrature)
      self.output.writeInfo()
      self.output.open(totalTime, numTimeSteps)

    self._logger.eventEnd(logEvent)
    return


  def poststep(self, t, dt, totalTime, fields):
    """
    Hook for doing stuff after advancing time step.
    """
    logEvent = "%spoststep" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    self._info.log("Writing fault data.")
    #self.output.writeData(t+dt, fields)

    self._logger.eventEnd(logEvent)
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


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    PetscComponent._configure(self)
    self.faultQuadrature = self.inventory.quadrature
    self.upDir = map(float, self.inventory.upDir)
    self.normalDir = map(float, self.inventory.normalDir)
    self.matDB = self.inventory.matDB
    ModuleFault.id(self, self.inventory.id)
    ModuleFault.label(self, self.inventory.label)
    return

  
  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    raise NotImplementedError("Please implement _createModuleObj() in " \
                              "derived class.")
  
  
# End of file 
