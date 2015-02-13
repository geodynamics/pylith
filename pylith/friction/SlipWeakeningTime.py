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

## @file pylith/friction/SlipWeakeningTime.py
##
## @brief Python object implementing slip weakening with forced
## weakening at a give time (sometimes used for nucleation).
##
## Factory: friction_model.

from FrictionModel import FrictionModel
from friction import SlipWeakeningTime as ModuleSlipWeakeningTime

# SlipWeakeningTime class
class SlipWeakeningTime(FrictionModel, ModuleSlipWeakeningTime):
  """
  Python object implementing Slip Weakening.

  Factory: friction_model.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="slipweakeningtime"):
    """
    Constructor.
    """
    FrictionModel.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': ["static_coefficient",
                     "dynamic_coefficient",
                     "slip_weakening_parameter",
                     "cohesion",
                     "weakening_time"],
            'data': ["cumulative_slip",
                     "previous_slip"]},
         'cell': \
           {'info': [],
            'data': []}}
    self._loggingPrefix = "FrSWT "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleSlipWeakeningTime.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def friction_model():
  """
  Factory associated with SlipWeakeningTime.
  """
  return SlipWeakeningTime()


# End of file 
