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

## @file pylith/materials/ElasticMaterial.py
##
## @brief Python abstract base class for managing physical properties
## of an elastic material.
##
## Factory: material

from Material import Material

# ElasticMaterial class
class ElasticMaterial(Material):
  """
  Python abstract base class for managing physical properties of an
  elastic material.

  Factory: material
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticmaterial"):
    """
    Constructor.
    """
    Material.__init__(self, name)
    return


  def useElasticBehavior(self, flag):
    """
    Set useElasticBehavior flag (True=elastic, False=inelastic if applicable).
    """
    assert(None != self.cppHandle)
    self.cppHandle.useElasticBehavior = flag
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Material._configure(self)
    return

  
# End of file 
