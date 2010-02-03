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

## @file pylith/friction/StaticFriction.py
##
## @brief Python object implementing static friction.
##
## Factory: friction_model.

from FrictionModel import FrictionModel
from friction import StaticFriction as ModuleStaticFriction

# StaticFriction class
class StaticFriction(FrictionModel, ModuleStaticFriction):
  """
  Python object implementing static friction.

  Factory: friction_model.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="staticfriction"):
    """
    Constructor.
    """
    FrictionModel.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': ["friction_coefficient"],
            'data': []},
         'cell': \
           {'info': [],
            'data': []}}
    self._loggingPrefix = "FrStat "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleStaticFriction.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def friction_model():
  """
  Factory associated with StaticFriction.
  """
  return StaticFriction()


# End of file 
