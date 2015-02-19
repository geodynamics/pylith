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

## @file pylith/friction/SlipWeakeningTimeStable.py
##
## @brief Python object implementing slip weakening with forced
## weakening at a give time (sometimes used for nucleation).
##
## Factory: friction_model.

from FrictionModel import FrictionModel
from friction import SlipWeakeningTimeStable as ModuleSlipWeakeningTimeStable

# SlipWeakeningTimeStable class
class SlipWeakeningTimeStable(FrictionModel, ModuleSlipWeakeningTimeStable):
  """
  Python object implementing Slip Weakening.

  Factory: friction_model.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="slipweakeningtimestable"):
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
                     "time_weakening_time",
                     "time_weakening_parameter"],
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
    ModuleSlipWeakeningTimeStable.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def friction_model():
  """
  Factory associated with SlipWeakeningTimeStable.
  """
  return SlipWeakeningTimeStable()


# End of file 
