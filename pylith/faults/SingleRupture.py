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

from pylith.utils.PetscComponent import PetscComponent

# SingleRupture class
class SingleRupture(PetscComponent):
  """
  Python kinematic rupture container with one rupture.

  Inventory

  @class Inventory
  Python object for managing SingleRupture facilities and properties.
  
  \b Properties
  @li None
  
  \b Facilities
  @li \b rupture Kinematic earthquake rupture in problem

  """

  # INVENTORY //////////////////////////////////////////////////////////

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
    PetscComponent.__init__(self, name, facility="rupture")
    return


# End of file 
