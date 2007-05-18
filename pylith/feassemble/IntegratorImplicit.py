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

## @file pylith/feassemble/IntegratorImplicit.py
##
## @brief Python object for implicit time integration of actions with
## finite-elements.
##
## Factory: integrator

from Integrator import Integrator

# IntegratorInertia class
class IntegratorImplicit(Integrator):
  """
  Python object for implicit integration of operator actions with
  finite-elements.

  Factory: integrator.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="integratorimplicit"):
    """
    Constructor.
    """
    Integrator.__init__(self, name)
    return


  def timeStep(self, t):
    """
    Set time step for advancing from time t to time t+dt.
    """
    self.cppHandle.timeStep = t.value
    return


  def stableTimeStep(self):
    """
    Get stable time step for advancing from time t to time t+dt.
    """
    return self.cppHandle.getStableTimeStep()


  def integrateResidual(self, fieldOut, fieldInT):
    """
    Integrate residual term for quasi-static terms for finite-elements.
    """
    self.cppHandle.integrateResidual(fieldOut, fieldInT, self.mesh.cppHandle)
    return


  def integrateJacobian(self, jacobian, fieldInT):
    """
    Integrate Jacobian term for quasi-static terms for finite-elements.
    """
    self.cppHandle.integrateJacobian(jacobian, fieldInT, self.mesh.cppHandle)
    return


# End of file 
