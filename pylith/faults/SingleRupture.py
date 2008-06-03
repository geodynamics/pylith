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
##
## Factory: object_bin

from pylith.utils.ObjectBin import ObjectBin

# SingleRupture class
class SingleRupture(ObjectBin):
  """
  Python kinematic rupture container with one rupture.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ObjectBin.Inventory):
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
    ObjectBin.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set attributes from inventory.
    """
    ObjectBin._configure(self)
    self.bin = [self.inventory.rupture]
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def object_bin():
  """
  Factory associated with Homogeneous.
  """
  return SingleRupture()


# End of file 
