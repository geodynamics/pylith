#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/friction/TimeWeakening.py
##
## @brief Python object implementing time Weakening.
##
## Factory: friction_model.

from FrictionModel import FrictionModel
from friction import TimeWeakening as ModuleTimeWeakening

# TimeWeakening class
class TimeWeakening(FrictionModel, ModuleTimeWeakening):
  """
  Python object implementing Time Weakening.

  Factory: friction_model.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timeweakening"):
    """
    Constructor.
    """
    FrictionModel.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': ["static_coefficient",
                     "dynamic_coefficient",
                     "time_weakening_parameter",
                     "cohesion"],
            'data': ["elapsed_time"]},
         'cell': \
           {'info': [],
            'data': []}}
    self._loggingPrefix = "FrTmWk "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleTimeWeakening.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def friction_model():
  """
  Factory associated with TimeWeakening.
  """
  return TimeWeakening()


# End of file 
