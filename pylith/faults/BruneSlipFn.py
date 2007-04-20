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

## @file pylith/faults/BruneSlipFn.py
##
## @brief Python object for slip time function that follows the
## integral of Brune's (1970) far-field time function.
##
## Factory: slip_time_fn

from SlipTimeFn import SlipTimeFn

# BruneSlipFn class
class BruneSlipFn(SlipTimeFn):
  """
  Python object for slip time function that follows the integral of
  Brune's (1970) far-field time function.

  Factory: slip_time_fn
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="bruneslipfn"):
    """
    Constructor.
    """
    SlipTimeFn.__init__(self, name)
    import pylith.faults.faults as bindings
    self.cppHandle = bindings.BruneSlipFn()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    SlipTimeFn._configure(self)
    return

  
# End of file 
