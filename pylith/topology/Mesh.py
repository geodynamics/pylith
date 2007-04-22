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

## @file pylith/topology/Mesh.py
##
## @brief Python Mesh for finite-element topology information.
##
## Factory: finite_element_mesh

from pyre.components.Component import Component

# Mesh class
class Mesh(Component):
  """
  Python Mesh for finite-element topology information.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="mesh"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="finite_element_mesh")
    import pylith.topology.topology as bindings
    self.cppHandle = bindings.Mesh()
    self.coordsys = None
    self.debug = False
    return


  def initialize(self, coordsys):
    """
    Initialize mesh.
    """
    self.coordsys = coordsys
    return


  def setDebug(self, flag):
    """
    Set debugging flag.
    """
    self.debug = flag
    self.cppHandle.debug = self.debug
    return
  

  def dimension(self):
    """
    Get dimension of mesh.
    """
    dim = None
    if not self.cppHandle is None:
      dim = self.cppHandle.dimension
    return dim

  def distribute(self):
    """
    Distribute mesh across processors.
    """
    self._info.log("WARNING: Mesh::distribute() not tested.")
    self.cppHandle.distribute()
    return self


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def finite_element_mesh():
  """
  Factory associated with Mesh.
  """
  return Mesh()


# End of file
