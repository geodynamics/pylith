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

## @file pylith/faults/EqKinSrc.py
##

## @brief Python object for managing parameters for a kinematic
## earthquake sources.
##
## EqKinSrc is responsible for providing the value of slip at time t
## over a fault surface.
##
## Factory: eq_kinematic_src

from Component import Component

# EqKinSrc class
class EqKinSrc(Component):
  """
  Python object for managing parameters for a kinematic earthquake sources.

  Factory: eq_kinematic_src
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing EqKinSrc facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing EqKinSrc facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b slip_function Analytical form for slip time function

    import pyre.inventory

    from BruneSlipFn import BruneSlipFn
    slipFn = pyre.inventory.facility("slip_function", family="slip_time_fn",
                                     factory=BruneSlipFn)
    slipFn.meta['tip'] = "Slip time history function."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="eqkinsrc"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="eqkinsrc")
    return


  def initialize(self, mesh):
    """
    Initialize cohesive elements.
    """
    Component.initialize(self, mesh)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    slipFn = self.inventory.slipFn
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def eq_kinematic_src():
  """
  Factory associated with EqKinSrc.
  """
  return EqKinSrc()


# End of file 
