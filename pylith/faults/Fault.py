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

from pyre.components.Component import Component

# Validator for up-direction
def validateUpDir(value):
  """
  Validate up-direction.
  """
  msg = "Up-direction must be a 3 component vector (list)."
  if not isinstance(value, list):
    raise ValueError(msg)
  if 3 != len(list):
    raise ValueError(msg)
  try:
    nums = map(float, value)
  except:
    raise ValueError(msg)
  return value


# Fault class
class Fault(Component):
  """
  Python abstract base class for a fault surface.

  This implementation of a fault associates both physical
  properties and a quadrature scheme with the fault.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Fault facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Fault facilities and properties.
    ##
    ## \b Properties
    ## @li \b id Fault identifier
    ## @li \b name Name of fault
    ## @li \b up-dir Up-dip or up direction
    ##   (perpendicular to along-strike and not collinear with fault normal)
    ##
    ## \b Facilities
    ## @li \b quadrature Quadrature object for numerical integration

    import pyre.inventory

    id = pyre.inventory.int("id", default=100)
    id.meta['tip'] = "Fault identifier (must be unique across all faults " \
                     "and materials)."

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Name of material."

    upDir = pyre.inventory.list("up-dir", default=[0, 0, 1],
                                validator=validateUpDir)
    upDir.meta['tip'] = "Up-dip or up direction " \
                        "(perpendicular to along-strike and not collinear " \
                        "with fault normal)."

    from pylith.feassemble.quadrature.Quadrature import Quadrature
    quadrature = pyre.inventory.facility("quadrature", factory=Quadrature)
    quadrature.meta['tip'] = "Quadrature object for numerical integration."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fault"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="fault")
    self.cppHandle = None
    return


  def adjustTopology(self, mesh):
    """
    Adjust mesh topology for fault implementation.
    """
    assert(None != self.cppHandle)
    self.cppHandle.id = self.id
    self.cppHandle.label = self.label
    self.cppHandle.adjustTopology(mesh.cppHandle)
    return


  def initialize(self, mesh):
    """
    Initialize fault.
    """
    self._info.log("Initializing fault '%s'." % self.label)

    self.quadrature.initialize()

    faultDim = mesh.dimension() - 1
    if faultDim != self.quadrature.cell.cellDim:
      raise ValueError, \
            "Quadrature is incompatible with fault surface.\n" \
            "Dimensions for quadrature: %d, dimensions of fault: %d" % \
            (self.quadrature.cell.cellDim, faultDim)

    assert(None != self.cppHandle)
    self.cppHandle.id = self.id
    self.cppHandle.label = self.label
    self.cppHandle.quadrature = self.quadrature.cppHandle
    self.cppHandle.initialize(mesh.cppHandle, mesh.coordsys.cppHandle,
                              self.upDir)
    return


  def timeStep(self, dt):
    """
    Set time step for advancing from time t to time t+dt.
    """
    assert(None != self.cppHandle)
    self.cppHandle.timeStep = dt.value
    return


  def stableTimeStep(self):
    """
    Get stable time step for advancing from time t to time t+dt.
    """
    assert(None != self.cppHandle)
    from pyre.units.time import second
    return self.cppHandle.stableTimeStep*second


  def integrateResidual(self, residual, fields):
    """
    Integrate contributions to residual term at time t.
    """
    assert(None != self.cppHandle)
    self.cppHandle.integrateResidual(residual, fields.cppHandle,
                                     self.mesh.cppHandle)
    return


  def needNewJacobian(self):
    """
    Returns true if we need to recompute Jacobian matrix for operator,
    false otherwise.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.needNewJacobian


  def integrateJacobian(self, jacobian, fields):
    """
    Integrate contributions to Jacobian term at time t.
    """
    assert(None != self.cppHandle)
    self.cppHandle.integrateJacobian(jacobian, fields.cppHandle,
                                     self.mesh.cppHandle)
    return


  def updateState(self, field):
    """
    Update state variables as needed.
    """
    assert(None != self.cppHandle)
    self.cppHandle.updateState(field, self.mesh.cppHandle)
    return
    

  def setConstraintSizes(self, field):
    """
    Set constraint sizes in field.
    """
    assert(None != self.cppHandle)
    self.cppHandle.setConstraintSizes(field, self.mesh.cppHandle)
    return


  def setConstraints(self, field):
    """
    Set constraints for field.
    """
    assert(None != self.cppHandle)
    self.cppHandle.setConstraints(field, self.mesh.cppHandle)
    return


  def setField(self, t, field):
    """
    Set constrained values in field at time t.
    """
    assert(None != self.cppHandle)
    self.cppHandle.setField(t.value, field, self.mesh.cppHandle)
    return


  def finalize(self):
    """
    Cleanup after time stepping.
    """
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    self.id = self.inventory.id
    self.label = self.inventory.label
    self.upDir = self.inventory.upDir
    self.quadrature = self.inventory.quadrature
    return

  
# End of file 
