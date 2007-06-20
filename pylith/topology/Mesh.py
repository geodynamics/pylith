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
    self.cppHandle = None
    self.coordsys = None
    self.debug = False
    return


  def initialize(self, coordsys):
    """
    Initialize mesh.
    """
    self._createCppHandle()
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
  

  def view(self):
    """
    View the mesh.
    """
    if not self.cppHandle is None:
      self.cppHandle.view
    return


  def comm(self):
    """
    Get MPI communicator associated with mesh.
    """
    comm = None
    #if not self.cppHandle is None:
    #  comm = self.cppHandle.comm
    import mpi
    comm = mpi.MPI_COMM_WORLD
    return comm


  def getRealSection(self, label):
    """
    Get real section from mesh.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.getRealSection(label)


  def createRealSection(self, label, fiberDim):
    """
    Create real section (but don't allocate).
    """
    assert(None != self.cppHandle)
    return self.cppHandle.createRealSection(label, fiberDim)


  def allocateRealSection(self, field):
    """
    Allocated (and zero) real section.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.allocateRealSection(field)


  def createMatrix(self, field):
    """
    Create sparse matrix compatible with field.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.createMatrix(field)
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    return

  
  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.topology.topology as bindings
      self.cppHandle = bindings.Mesh()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def finite_element_mesh():
  """
  Factory associated with Mesh.
  """
  return Mesh()


# End of file
