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

## @file pylith/problems/SolnDispVel.py
##
## @brief Python subfields container with displacement and velocity subfields.

from pylith.utils.PetscComponent import PetscComponent

# SolnDispVel class
class SolnDispVel(PetscComponent):
  """
  Python subfields container with displacement and velocity subfields.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing SolnDispVel facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Homogeneous facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b displacement Displacement subfield.
    ## @li \b velocity Velicity subfield.

    import pyre.inventory

    from SubfieldDisplacement import SubfieldDispacement
    displacement = pyre.inventory.facility("displacement", family="subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement field."

    from SubfieldVelocity import SubfieldVelocity
    velocity = pyre.inventory.facility("velocity", family="subfield", factory=SubfieldVelocity)
    velocity.meta['tip'] = "Velocity field."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solndispvel"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="subfields")
    return


# End of file 
