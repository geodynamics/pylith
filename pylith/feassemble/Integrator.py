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

## @file pylith/feassemble/Integrator.py
##
## @brief Python abstract base class for integration of operator
## actions with finite-elements.
##
## Factory: fe_integrator.

def implementsIntegrator(obj):
  """
  Check whether object implements an integrator.
  """
  result = True
  attrs = dir(obj)
  if not "timeStep" in attrs or \
     not "stableTimeStep" in attrs or \
     not "useSolnIncr" in attrs or \
     not "integrateResidual" in attrs or \
     not "integrateJacobian" in attrs or \
     not "updateState" in attrs or \
     not "verifyConfiguration" in attrs or \
     not "finalize" in attrs:
    result = False
  return result


# Integrator class
class Integrator(object):
  """
  Python abstract base class for integration of actions with
  finite-elements.

  Factory: integrator.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    self.quadrature = None
    self.mesh = None
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    assert(None != self.cppHandle)
    self.cppHandle.verifyConfiguration(self.mesh.cppHandle)
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


  def useSolnIncr(self, flag):
    """
    Set behavior for using total field solution or incremental field solution.
    """
    assert(None != self.cppHandle)
    self.cppHandle.useSolnIncr = flag
    return
  

  def integrateResidual(self, residual, t, fields):
    """
    Integrate contributions to residual term at time t.
    """
    assert(None != self.cppHandle)
    self.cppHandle.integrateResidual(residual, t.value, fields.cppHandle,
                                     self.mesh.cppHandle)
    return


  def needNewJacobian(self):
    """
    Returns true if we need to recompute Jacobian matrix for operator,
    false otherwise.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.needNewJacobian


  def integrateJacobian(self, jacobian, t, fields):
    """
    Integrate contributions to Jacobian term at time t.
    """
    assert(None != self.cppHandle)
    self.cppHandle.integrateJacobian(jacobian, t.value, fields.cppHandle,
                                     self.mesh.cppHandle)
    return


  def updateState(self, t, field):
    """
    Update state variables as needed.
    """
    assert(None != self.cppHandle)
    self.cppHandle.updateState(t.value, field, self.mesh.cppHandle)
    return
    

  def finalize(self):
    """
    Cleanup after time stepping.
    """
    return
  

# End of file 
