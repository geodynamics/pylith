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


  def integrateConstant(self, fieldOut, fieldInT, fieldInTmdt):
    """
    Integrate constant term for dynamic terms for finite-elements.
    """
    self.cppHandle.integrateConstant(fieldOut, fieldInT, fieldInTmdt,
                                     self.mesh.cppHandle)
    return


  def integrateJacobian(self, jacobian, fieldInT, coords):
    """
    Integrate Jacobian term for dynamic terms for finite-elements.
    """
    self.cppHandle.integrateJacobian(jacobian, fieldInT, self.mesh.cppHandle)
    return


# End of file 
