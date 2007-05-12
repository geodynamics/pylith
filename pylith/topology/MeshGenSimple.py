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

## @file pylith/topology/MeshGenSimple.py
##
## @brief Python manager for simple mesh generator.
##
## Factory: mesh_generator.

from MeshGenerator import MeshGenerator

# MeshGenSimple class
class MeshGenSimple(MeshGenerator):
  """
  Python manager for simple mesh generator.

  Factory: mesh_generator
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshgensimple"):
    """
    Constructor.
    """
    MeshGenerator.__init__(self, name)
    import pylith.topology.topology as bindings
    self.cppHandle = bindings.MeshGenSimple()
    return


  def create(self, faults=None):
    """
    Generate a Mesh from a boundary
    """
    return self.cppHandle.generate(self.boundary)


  def setBoundary(self, boundary):
    """
    Set boundary for domain to mesh.
    """
    self.boundary = boundary
    return
  

  def createCubeBoundary(self):
    """
    Returns a Mesh that is the boundary of the unit cube
    """
    return self.cppHandle.createCubeBoundary(self.debug)


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    MeshGenerator._configure(self)
    return


# End of file 
