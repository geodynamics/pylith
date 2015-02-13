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

## @file pylith/meshio/SingleOutput.py
##
## @brief Python container with one output manager.

from pylith.utils.PetscComponent import PetscComponent

# SingleOutput class
class SingleOutput(PetscComponent):
  """
  Python container with one output manager.

  Factory: object_bin
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing SingleOutput facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing SingleOutput facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b output Output manager

    import pyre.inventory

    from OutputSoln import OutputSoln
    output = pyre.inventory.facility("output", family="output_manager",
                                     factory=OutputSoln)
    output.meta['tip'] = "Output manager."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="singleoutput"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="output")
    return


# End of file 
