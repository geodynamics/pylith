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

## @file pylith/faults/SingleRupure.py
##
## @brief Python kinematic rupture container with one rupture.

from pyre.components.Component import Component

# SingleRupture class
class SingleRupture(Component):
  """
  Python kinematic rupture container with one rupture.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing SingleRupture facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing SingleRupture facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b rupture Kinematic earthquake rupture in problem

    import pyre.inventory

    from EqKinSrc import EqKinSrc
    rupture = pyre.inventory.facility("rupture", family="eq_kinematic_src",
                                       factory=EqKinSrc)
    rupture.meta['tip'] = "Kinematic earthquake rupture in problem."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="singlerupture"):
    """
    Constructor.
    """
    Component.__init__(self, name)
    return


# End of file 
