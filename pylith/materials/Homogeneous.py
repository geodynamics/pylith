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

## @file pylith/materials/Homogeneous.py
##
## @brief Python materials container with one material.

from pylith.utils.PetscComponent import PetscComponent

# Homogeneous class
class Homogeneous(PetscComponent):
  """
  Python materials container with one material.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing Homogeneous facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Homogeneous facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b material Material in problem

    import pyre.inventory

    from ElasticIsotropic3D import ElasticIsotropic3D
    material = pyre.inventory.facility("material", family="material",
                                       factory=ElasticIsotropic3D)
    material.meta['tip'] = "Material in problem."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="homogeneous"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="material")
    return


# End of file 
