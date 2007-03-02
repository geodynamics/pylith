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

## @file pylith/topology/MeshGenerator.py
##
## @brief Python abstract base class for mesh generator.
##
## Factory: mesh_generator.

from pyre.components.Component import Component

# MeshGenerator class
class MeshGenerator(Component):
  """
  Python abstract base class for mesh generator.

  Factory: mesh_generator
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshgenerator"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="mesh_generator")
    return


  def create(self):
    """
    Hook for creating mesh.
    """
    raise NotImplementedError, "MeshGenerator::create() not implemented."
    return


# End of file 
