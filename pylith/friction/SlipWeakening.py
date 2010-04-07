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

## @file pylith/friction/SlipWeakening.py
##
## @brief Python object implementing slip Weakening.
##
## Factory: friction_model.

from FrictionModel import FrictionModel
from friction import SlipWeakening as ModuleSlipWeakening

# SlipWeakening class
class SlipWeakening(FrictionModel, ModuleSlipWeakening):
  """
  Python object implementing Slip Weakening.

  Factory: friction_model.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="slipweakening"):
    """
    Constructor.
    """
    FrictionModel.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': ["static_coefficient","dynamic_coefficient","slip_weakening_parameter","cohesion"],
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
    ModuleSlipWeakening.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def friction_model():
  """
  Factory associated with SlipWeakening.
  """
  return SlipWeakening()


# End of file 
