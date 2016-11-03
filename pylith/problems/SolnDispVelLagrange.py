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
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/problems/SolnDispVelLagrange.py
##
## @brief Python subfields container with displacement, velocity, and
## fault Lagrange multiplier subfields.

from pylith.utils.PetscComponent import PetscComponent

# SolnDispVelLagrange class
class SolnDispVelLagrange(PetscComponent):
  """
  Python subfields container with displacement, velocity, and fault
  Lagrange multiplier subfields.

  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing SolnDispVelLagrange facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing SolnDispVelLagrange facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b displacement Displacement subfield.
    ## @li \b velocity Velocity subfield.
    ## @li \b lagrange_fault Fault Lagrange multiplier subfield.

    import pyre.inventory

    from SubfieldDisplacement import SubfieldDisplacement
    displacement = pyre.inventory.facility("displacement", family="subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    from SubfieldVelocity import SubfieldVelocity
    velocity = pyre.inventory.facility("velocity", family="subfield", factory=SubfieldVelocity)
    velocity.meta['tip'] = "Velocity subfield."

    from SubfieldLagrangeFault import SubfieldLagrangeFault
    lagrangeFault = pyre.inventory.facility("lagrange_fault", family="subfield", factory=SubfieldLagrangeFault)
    lagrangeFault.meta['tip'] = "Fault Lagrange multiplier subfield."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solndispvel"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="subfields")
    return


  def _configure(self):
    PetscComponent._configure(self)
    self.displacement = self.inventory.displacement
    self.velocity = self.inventory.velocity
    self.lagrangeFault = self.inventory.lagrangeFault
    return


  def components(self):
    """
    Order of facilities in Inventory is ambiguous, so overwrite
    components() to insure order is [displacement, velocity, lagrange_fault].

    """
    return [self.inventory.displacement, self.inventory.velocity, self.inventory.lagrangeFault]


# End of file 
