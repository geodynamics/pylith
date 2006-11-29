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

## @file pylith/feassemble/IntegratorInertia.py

## @brief Python object for integration of inertial operator
## actions with finite-elements.

from Integrator import Integrator

# IntegratorInertia class
class IntegratorInertia(Integrator):
  """
  Python object for integration of inertial operator actions with
  finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="integratorinertia"):
    """
    Constructor.
    """
    Integrator.__init__(self, name)

    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.IntegratorInertia()
    return


# End of file 
