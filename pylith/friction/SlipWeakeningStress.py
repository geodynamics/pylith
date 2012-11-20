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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/friction/SlipWeakeningStress.py
##
## @brief Python object implementing slip weakening with stress levels
## and forced weakening at a give time (sometimes used for
## nucleation).
##
## Factory: friction_model.

from FrictionModel import FrictionModel
from friction import SlipWeakeningStress as ModuleSlipWeakeningStress

# SlipWeakeningStress class
class SlipWeakeningStress(FrictionModel, ModuleSlipWeakeningStress):
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
           {'info': ["static_stress",
                     "dynamic_stress",
                     "slip_weakening_parameter",
                     "weakening_time"],
            'data': ["cumulative_slip",
                     "previous_slip"]},
         'cell': \
           {'info': [],
            'data': []}}
    self._loggingPrefix = "FrSWS "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleSlipWeakeningStress.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def friction_model():
  """
  Factory associated with SlipWeakeningStress.
  """
  return SlipWeakeningStress()


# End of file 
