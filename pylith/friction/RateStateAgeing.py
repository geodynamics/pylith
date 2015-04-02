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

## @file pylith/friction/RateStateAgeing.py
##
## @brief Python object implementing Rate and State with Ageing Law.
##
## Factory: friction_model.

from FrictionModel import FrictionModel
from friction import RateStateAgeing as ModuleRateStateAgeing

# RateStateAgeing class
class RateStateAgeing(FrictionModel, ModuleRateStateAgeing):
  """
  Python object implementing Rate and State with Ageing Law.

  Factory: friction_model.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(FrictionModel.Inventory):
    """
    Python object for managing RateStateAgeing facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing RateStateAgeing facilities and properties.
    ##
    ## \b Properties
    ## @li \b linear_slip_rate Nondimensional slip rate below which friction 
    ## varies linearly with slip rate.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    linearSlipRate = pyre.inventory.float("linear_slip_rate", default=1.0e-12,
                                          validator=pyre.inventory.greaterEqual(0.0))
    linearSlipRate.meta['tip'] = "Nondimensional slip rate below which friction " \
        "varies linearly with slip rate."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="ratestateageing"):
    """
    Constructor.
    """
    FrictionModel.__init__(self, name)
    self.availableFields = \
        {'vertex': \
           {'info': ["reference_friction_coefficient",
                     "reference_slip_rate",
                     "characteristic_slip_distance",
                     "constitutive_parameter_a",
                     "constitutive_parameter_b",
                     "cohesion"],
            'data': ["state_variable"]},
         'cell': \
           {'info': [],
            'data': []}}
    self._loggingPrefix = "FrRSAg "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    try:
      FrictionModel._configure(self)
      ModuleRateStateAgeing.linearSlipRate(self, self.inventory.linearSlipRate)
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring friction model "
                       "(%s):\n%s" % (aliases, err.message))
    return

  
  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    ModuleRateStateAgeing.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def friction_model():
  """
  Factory associated with RateStateAgeing.
  """
  return RateStateAgeing()


# End of file 
