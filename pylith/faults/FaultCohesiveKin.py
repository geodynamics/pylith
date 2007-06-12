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

## @file pylith/faults/FaultCohesiveKin.py
##

## @brief Python object for a fault surface with kinematic
## (prescribed) slip implemented with cohesive elements.
##
## Factory: fault

from FaultCohesive import FaultCohesive

# FaultCohesiveKin class
class FaultCohesiveKin(FaultCohesive):
  """
  Python object for a fault surface with kinematic (prescribed) slip
  implemented with cohesive elements.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(FaultCohesive.Inventory):
    """
    Python object for managing FaultCohesiveKin facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing FaultCohesiveKin facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b eq_src Kinematic earthquake source information.

    import pyre.inventory

    from EqKinSrc import EqKinSrc
    eqsrc = pyre.inventory.facility("eq_src", family="eq_kinematic_src",
                                    factory=EqKinSrc)
    eqsrc.meta['tip'] = "Kinematic earthquake source information."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivekin"):
    """
    Constructor.
    """
    FaultCohesive.__init__(self, name)
    import pylith.faults.faults as bindings
    self.cppHandle = bindings.FaultCohesiveKin()
    return


  def initialize(self, mesh):
    """
    Initialize cohesive elements.
    """
    self.mesh = mesh
    assert(None != self.cppHandle)
    self.eqsrc.initialize()
    self.cppHandle.eqsrc = self.eqsrc.cppHandle
    FaultCohesive.initialize(self, mesh)
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
    FaultCohesive._configure(self)
    self.eqsrc = self.inventory.eqsrc
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def fault():
  """
  Factory associated with FaultCohesiveKin.
  """
  return FaultCohesiveKin()


# End of file 
