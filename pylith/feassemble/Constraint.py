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

## @file pylith/feassemble/Constraint.py
##
## @brief Python abstract base class for constraints on operator
## actions with finite-elements.
##
## Factory: fe_constraint.

def implementsConstraint(obj):
  """
  Check whether object implements a constraint.
  """
  result = True
  attrs = dir(obj)
  if not "timeStep" in attrs or \
     not "setConstraintSizes" in attrs or \
     not "setConstraints" in attrs or \
     not "setField" in attrs or \
     not "finalize" in attrs:
    result = False
  return result


# Constraint class
class Constraint(object):
  """
  Python abstract base class for constraints on operator
  actions with finite-elements.

  Factory: constraint.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self):
    """
    Constructor.
    """
    self.mesh = None
    return


  def timeStep(self, dt):
    """
    Set time step for advancing from time t to time t+dt.
    """
    assert(None != self.cppHandle)
    self.cppHandle.timeStep = dt.value
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
  

# End of file 
