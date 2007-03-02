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

## @file pylith/feassemble/IntegratorExplicit.py
##
## @brief Python object for explicit time integration of actions with
## finite-elements.
##
## Factory: integrator

from Integrator import Integrator

# IntegratorInertia class
class IntegratorExplicit(Integrator):
  """
  Python object for explicit integration of operator actions with
  finite-elements.

  Factory: integrator.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="integratorexplicit"):
    """
    Constructor.
    """
    Integrator.__init__(self, name)
    return


  def timeStep(self, t):
    """
    Set time step for advancing from time t to time t+dt.
    """
    self.cppHandle.timeStep(t.value)
    return


  def stableTimeStep(self):
    """
    Get stable time step for advancing from time t to time t+dt.
    """
    return self.cppHandle.getStableTimeStep()


  def integrateResidual(self,
                        residual, fieldInT, fieldInTmdt, coords, lumpJacobian):
    """
    Integrate residual term for dynamic terms for finite-elements.
    """
    if lumpJacobian:
      self.cppHandle.integrateResidualLumped(residual,
                                             fieldInT, fieldInTmdt, coords)
    else:
      self.cppHandle.integrateResidual(residual, fieldInT, fieldInTmdt, coords)
    return


  def integrateJacobian(self, jacobian, fieldInT, coords, lumpJacobian):
    """
    Integrate Jacobian term for dynamic terms for finite-elements.
    """
    if lumpJacobian:
      self.cppHandle.integrateJacobianLumped(jacobian, fieldInT, coords)
    else:
      self.cppHandle.integrateJacobian(jacobian, fieldInT, coords)
    return


# End of file 
