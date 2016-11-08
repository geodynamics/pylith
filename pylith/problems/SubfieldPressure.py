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

## @file pylith/problems/SubfieldPressure.py
##
## @brief Python object for pressure subfield.
##
## Factory: subfield.

from .SolutionSubfield import SolutionSubfield

# SubfieldPressure class
class SubfieldPressure(SolutionSubfield):
  """
  Python object for pressure subfield.

  Factory: subfield.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(SolutionSubfield.Inventory):
      """
      Python object for managing SubfieldPressure facilities and properties.
      """

    ## @class Inventory
    ## Python object for managing SubfieldPressure facilities and properties.
    ##
    ## \b Properties
    ## @li \b name Name for subfield.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    from .SolutionSubfield import validateName
    name = pyre.inventory.str("name", default="pressure", validator=validateName)
    name.meta['tip'] = "Name for subfield."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="subfieldpressure"):
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
    self.vectorFieldType = Field.SCALAR
    self.ncomponents = 1
    self.scale = normalizer.pressureScale()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def soln_subfield():
  """
  Factory associated with SubfieldPressure.
  """
  return SubfieldPressure()


# End of file
