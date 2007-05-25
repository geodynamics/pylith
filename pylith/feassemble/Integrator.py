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
  if not "initQuadrature" in attrs or \
     not "timeStep" in attrs or \
     not "stableTimeStep" in attrs or \
     not "integrateResidual" in attrs or \
     not "integrateJacobian" in attrs:
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


  def setMesh(self, mesh):
    """
    Set mesh.
    """
    self.mesh = mesh
    return
  

  def initQuadrature(self, quadrature):
    """
    Initialize quadrature.
    """
    assert(None != self.cppHandle)
    quadrature.initialize()
    self.quadrature = quadrature
    self.cppHandle.quadrature = self.quadrature.cppHandle
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
    return self.cppHandle.getStableTimeStep()


  def integrateResidual(self, residual, fields):
    """
    Integrate contributions to residual term at time t.
    """
    assert(None != self.cppHandle)
    self.cppHandle.integrateResidual(residual, fields, self.mesh.cppHandle)
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
    self.cppHandle.integrateJacobian(jacobian, fields, self.mesh.cppHandle)
    return


  def updateState(self, field):
    """
    Update state variables as needed.
    """
    assert(None != self.cppHandle)
    self.cppHandle.updateState(field, self.mesh.cppHandle)
    return
    

# End of file 
