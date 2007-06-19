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
from Mesh import Mesh

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
    self.cppHandle = None
    return


  def create(self, faults=None):
    """
    Generate a Mesh from a boundary
    """
    mesh = Mesh()
    self._createCppHandle()
    mesh.cppHandle = self.cppHandle.generate(self.boundary.cppHandle)
    return mesh


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
    mesh = Mesh()
    self._createCppHandle()
    mesh.cppHandle = self.cppHandle.createCubeBoundary(self.debug)
    return mesh


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    MeshGenerator._configure(self)
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.topology.topology as bindings
      self.cppHandle = bindings.MeshGenSimple()
    return
  

# End of file 
