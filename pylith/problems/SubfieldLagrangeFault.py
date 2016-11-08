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

## @file pylith/problems/SubfieldLagrangeFault.py
##
## @brief Python object for fault Lagrange multipliers subfield.
##
## Factory: subfield.

from .SolutionSubfield import SolutionSubfield

# SubfieldLagrangeFault class
class SubfieldLagrangeFault(SolutionSubfield):
  """
  Python object for fault Lagrange multipliers subfield.

  Factory: subfield.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(SolutionSubfield.Inventory):
    """
    Python object for managing SubfieldLagrangeFault facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing SubfieldLagrangeFault facilities and properties.
    ##
    ## \b Properties
    ## @li \b name Name for subfield.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    from .SolutionSubfield import validateName
    name = pyre.inventory.str("name", default="lagrange_multiplier_fault", validator=validateName)
    name.meta['tip'] = "Name for subfield."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="subfieldlagrangefault"):
    """
    Constructor.
    """
    SolutionSubfield.__init__(self, name)
    return


  def initialize(self, normalizer, spaceDim):
    """
    Initialize subfield metadata.
    """
    from pylith.topology.Field import Field
    self.vectorFieldType = Field.VECTOR
    self.ncomponents = spaceDim
    self.scale = normalizer.pressureScale()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def soln_subfield():
  """
  Factory associated with SubfieldLagrangeFault.
  """
  return SubfieldLagrangeFault()


# End of file
