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

## @file pylith/problems/EqDeformation.py
##
## @brief Python EqDeformation for computing deformation associated
## with earthquakes.
##
## Factory: problem.

from TimeDependent import TimeDependent

# EqDeformation class
class EqDeformation(TimeDependent):
  """
  Python EqDeformation for computing deformation associated with
  earthquakes.

  Factory: problem.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(TimeDependent.Inventory):
    """
    Python object for managing EqDeformation facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing EqDeformation facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b faults Faults or interior slip surfaces.

    import pyre.inventory

    from pylith.faults.FaultsBin import FaultsBin
    faults = pyre.inventory.facility("faults", family="faults",
                                     factory=FaultsBin)
    faults.meta['tip'] = "Faults or interior slip surfaces."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="eqdeformation"):
    """
    Constructor.
    """
    TimeDependent.__init__(self, name)
    return


  def checkpoint(self):
    """
    Save problem state for restart.
    """
    TimeDependent.checkpoint() # Save state of parent
    
    # Save state of this object
    raise NotImplementedError, \
          "EqDeformation::checkpoint() not implemented."

    # Save state of children
    self.faults.checkpoint()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    TimeDependent._configure(self)
    self.faults = self.inventory.faults
    return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
  """
  Factory associated with EqDeformation.
  """
  return EqDeformation()


# End of file 
