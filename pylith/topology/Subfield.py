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

## @file pylith/topology/Subfield.py
##
## @brief Python object for defining attributes of a subfield within a
## field.
##
## Factory: subfield.

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


# Subfield class
class Subfield(PetscComponent):
  """
  Python object for defining attributes of a subfield within a field.

  Factory: subfield.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing Subfield facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Subfield facilities and properties.
    ##
    ## \b Properties
    ## @li \b name Name for subfield.
    ## @li \b basis_order Order of basis functions.
    ## @li \b quad_order Order of numerical quadrature.
    ## @li \b basis_continuous Is basis continuous?
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory
    
    name = pyre.inventory.str("name", default="", validator=validateName)
    name.meta['tip'] = "Name for subfield."

    basisOrder = pyre.inventory.int("basis_order", default=1)
    basisOrder.meta['tip'] = "Order of basis functions."

    quadOrder = pyre.inventory.int("quad_order", default=1)
    quadOrder.meta['tip'] = "Order of numerical quadrature."

    basisContinuous = pyre.inventory.bool("basis_continous", default=True)
    basisContinuous.meta['tip'] = "Is basis continuous?"


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="subfield"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="subfield")

    # Set in derived class initialize().
    self.ncomponents = None
    self.vectorFieldType = None
    self.scale = None
    return


  def initialize(self, normalizer, spaceDim):
    """
    Initialize subfield metadata.
    """
    raise NotImplementedError("Implement in derived class.")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.name = self.inventory.name
    self.basisOrder = self.inventory.basisOrder
    self.quadOrder = self.inventory.quadOrder
    self.basisContinous = self.inventory.basisContinuous
    return


# FACTORIES ////////////////////////////////////////////////////////////

def subfield():
  """
  Factory associated with Subfield.
  """
  return Subfield()


# End of file 
