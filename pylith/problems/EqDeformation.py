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

## @brief Python EqDeformation for computing deformation associated
## with earthquakes.

from Problem import Problem

# EqDeformation class
class EqDeformation(Problem):
  """
  Python EqDeformation for computing deformation associated with
  earthquakes.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Problem.Inventory):
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

    #from Faults import Faults
    #faults = pyre.inventory.facility("faults", factory=Faults)
    #faults.meta['tip'] = "Faults or interior slip surfaces."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    """
    Setup formulation for solving PDE.
    """
    Problem.initialize(self)

    self._info.log("WARNING: EqDeofmraiton::initialize not implemented.")
    return


  def prestep(self):
    """
    User hook for doing stuff before advancing time step.
    """
    self._info.log("WARNING: EqDeformation::prestep not implemented.")
    return


  def step(self, dt):
    """
    Advance to next time step.
    """
    self._info.log("WARNING: EqDeformation::step not implemented.")
    return


  def poststep(self):
    """
    Update time and storage.
    """
    self._info.log("WARNING: EqDeformation::poststep not implemented.")
    return


  def stableTimestep(self):
    """
    Determine stable time step for problem.
    """
    self._info.log("WARNING: EqDeformation::stableTimestep not implemented.")
    return


  def checkpoint(self):
    """
    Save problem state for restart.
    """
    self._info.log("WARNING: EqDeformation::checkpoint not implemented.")
    return
  

  def __init__(self, name="eqdeformation"):
    """
    Constructor.
    """
    Problem.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Problem._configure(self)
    #self.faults = self.inventory.faults
    return


# End of file 
