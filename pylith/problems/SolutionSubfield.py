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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/problems/SolutionSubfield.py
##
## @brief Python solution field for problem.
##
## Factory: solution.

from pylith.utils.PetscComponent import PetscComponent

# Validator for name
def validateName(value):
  """
  Validate name of subfield.
  """
  if 0 == len(value):
    raise ValueError("Name of subfield not specified.")
  import re
  if re.search(r"\s", value):
    raise ValueError("Name of subfield cannot contain whitespace.")
  return value


# SolutionSubfield class
class SolutionSubfield(PetscComponent):
  """
  Python subfield in solution.

  Factory: solution.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing SolutionSubfield facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing SolutionSubfield facilities and properties.
    ##
    ## \b Properties
    ## @li \b components Number of components.
    ## @li \b vector_field_type Type of vector field ['scalar','vector','tensor'].
    ## @li \b scale Nondimensional scale for field.
    ## @li \b basis_order Order of basis functions.
    ## @li \b quad_order Order of numerical quadrature.
    ## @li \b basis_continuous Is basis continuous?
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory
    
    name = pyre.inventory.str("name", default="", validator=validateName)
    name.meta['tip'] = "Name for subfield."

    components = pyre.inventory.int("components", default=3)
    components.meta['tip'] = "Number of components."

    vectorFieldType = pyre.inventory.str("vector_field_type", default="vector", validator=pyre.inventory.choice(["scalar","vector","tensor"]))
    vectorFieldType.meta['tip'] = "Type of vector field ['scalar','vector','tensor']."

    # :TODO: Scale for field. Use Pyre units to make this scale in terms of nondimensional scales (pressure, length, time, etc).

    basisOrder = pyre.inventory.int("basis_order", default=1)
    basisOrder.meta['tip'] = "Order of basis functions."

    quadOrder = pyre.inventory.int("quad_order", default=1)
    quadOrder.meta['tip'] = "Order of numerical quadrature."

    basisContinuous = pyre.inventory.bool("basis_continous", default=True)
    basisContinuous.meta['tip'] = "Is basis continuous?"


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solution"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="solution")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.name = self.inventory.name
    self.components = self.inventory.components
    self.basisOrder = self.inventory.basisOrder
    self.quadOrder = self.inventory.quadOrder
    self.basisContinous = self.inventory.basisContinuous
    return


# FACTORIES ////////////////////////////////////////////////////////////

def subfield():
  """
  Factory associated with SolutionSubfield.
  """
  return SolutionSubfield()


# End of file 
