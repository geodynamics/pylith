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

## @file pylith/faults/FaultCohesiveKin.py
##

## @brief Python object for a fault surface with kinematic
## (prescribed) slip implemented with cohesive elements.
##
## Factory: fault

from FaultCohesive import FaultCohesive

# FaultCohesiveKin class
class FaultCohesiveKin(FaultCohesive):
  """
  Python object for a fault surface with kinematic (prescribed) slip
  implemented with cohesive elements.

  Factory: fault
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(FaultCohesive.Inventory):
    """
    Python object for managing FaultCohesiveKin facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing FaultCohesiveKin facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b eq_src Kinematic earthquake source information.

    import pyre.inventory

    from EqKinSrc import EqKinSrc
    eqSrc = pyre.inventory.facility("eq_src", family="eq_kinematic_src",
                                    factory=EqKinSrc)
    eqSrc.meta['tip'] = "Kinematic earthquake source information."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="faultcohesivekin"):
    """
    Constructor.
    """
    FaultCohesive.__init__(self, name)
    import pylith.faults.faults as bindings
    self.cppHandle = bindings.FaultCohesiveKin()
    return


  def initialize(self, mesh):
    """
    Initialize cohesive elements.
    """
    FaultCohesive.initialize(self, mesh)
    self.slip.initialize()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    FaultCohesive._configure(self)
    eqSrc = self.inventory.eqSrc
    return

  
# FACTORIES ////////////////////////////////////////////////////////////

def fault():
  """
  Factory associated with FaultCohesiveKin.
  """
  return FaultCohesiveKin()


# End of file 
