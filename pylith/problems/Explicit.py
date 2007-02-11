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

## @file pylith/problems/Explicit.py

## @brief Python Explicit object for solving equations using an
## explicit formulation.

from Formulation import Formulation

# Explicit class
class Explicit(Formulation):
  """
  Python Explicit object for solving equations using an explicit
  formulation.

  The formulation has the general form, [A(t)] {u(t+dt)} = {b(t)},
  where we want to solve for {u(t+dt)}, A(t) is usually constant
  (i.e., independent of time), and {b(t)} usually depends on {u(t)}
  and {u(t-dt)}.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Formulation.Inventory):
    """
    Python object for managing Explicit facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Explicit facilities and properties.
    ##
    ## \b Properties
    ## @li \b lump_jacobian Flag for indicating to use lumped
    ##   formulation for Jacobian matrix
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    lumpJacobian = pyre.inventory.bool("lump_jacobian", default=True)
    lumpJacobian.meta['tip'] = "Flag for indicating to use lumped " \
                               "formulation for Jacobian matrix."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def stableTimeStep(self):
    """
    Get stable time step for advancing forward in time.
    """
    self._info.log("WARNING: Explicit::stableTimeStep() not implemented.")
    from pyre.units.time import second
    dt = 0.0*second
    return dt
  

  def initialize(self):
    """
    Initialize explicit formulation.
    """
    Formulation.initialize(self)
    return


  def calcResidual(self):
    """
    Compute residual, {b(t)}.
    """
    self._info.log("WARNING: Explicit::calcResidual not implemented.")
    return


  def calcJacobian(self):
    """
    Compute Jacobian, [A(t)].
    """
    self._info.log("WARNING: Explicit::calcJacobian not implemented.")
    return


  def __init__(self, name="explicit"):
    """
    Constructor.
    """
    Formulation.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Formulation._configure(self)
    self.lumpJacobian = self.inventory.lumpJacobian
    return


# End of file 
