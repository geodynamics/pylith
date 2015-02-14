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

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(FrictionModel.Inventory):
    """
    Python object for managing SlipWeakening facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing FrictionModel facilities and properties.
    ##
    ## \b Properties
    ## @li \b force_healing Force healing after every time step.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    forceHealing = pyre.inventory.bool("force_healing", default=False)
    forceHealing.meta['tip'] = "Force healing after every time step."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="slipweakening"):
    """
    Constructor.
    """
    FrictionModel.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': ["static_coefficient",
                     "dynamic_coefficient",
                     "slip_weakening_parameter",
                     "cohesion"],
            'data': ["cumulative_slip",
                     "previous_slip"]},
         'cell': \
           {'info': [],
            'data': []}}
    self._loggingPrefix = "FrSlWk "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    try:
      FrictionModel._configure(self)
      ModuleSlipWeakening.forceHealing(self, self.inventory.forceHealing)
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring friction model "
                       "(%s):\n%s" % (aliases, err.message))
    return

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
