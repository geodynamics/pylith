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

## @file pylith/problems/SolnDispPres.py
##
## @brief Python subfields container with displacement and pressure subfields.

from pylith.utils.PetscComponent import PetscComponent

# SolnDispPres class
class SolnDispPres(PetscComponent):
  """
  Python subfields container with displacement and pressure subfields.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing SolnDispPres facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Homogeneous facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b displacement Displacement subfield.
    ## @li \b pressure Pressure subfield.

    import pyre.inventory

    from SubfieldDisplacement import SubfieldDisplacement
    displacement = pyre.inventory.facility("displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    from SubfieldPressure import SubfieldPressure
    pressure = pyre.inventory.facility("pressure", family="soln_subfield", factory=SubfieldPressure)
    pressure.meta['tip'] = "Pressure subfield."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solndisppres"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="soln_subfields")
    return


  def _configure(self):
    PetscComponent._configure(self)
    self.displacement = self.inventory.displacement
    self.pressure = self.inventory.pressure
    return


  def components(self):
    """
    Order of facilities in Inventory is ambiguous, so overwrite
    components() to insure order is [displacement, pressure].

    """
    return [self.inventory.displacement, self.inventory.pressure]


# End of file 
