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

## @file pylith/friction/ViscousFriction.py
##
## @brief Python object implementing viscous friction.
##
## Factory: friction_model.

from FrictionModel import FrictionModel
from friction import ViscousFriction as ModuleViscousFriction

# ViscousFriction class
class ViscousFriction(FrictionModel, ModuleViscousFriction):
  """
  Python object implementing viscous friction.

  Factory: friction_model.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="viscousfriction"):
    """
    Constructor.
    """
    FrictionModel.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': ["static_coefficient",
                     "reference_slip_rate"],
            'data': ["slip_rate"]},
         'cell': \
           {'info': [],
            'data': []}}
    self._loggingPrefix = "FrVisc "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleViscousFriction.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def friction_model():
  """
  Factory associated with ViscousFriction.
  """
  return ViscousFriction()


# End of file 
