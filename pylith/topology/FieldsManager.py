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

## @file pylith/topology/FieldsManager.py
## 
## @brief Python manager for fields over a mesh.

# FieldsManager class
class FieldsManager(object):
  """
  Python manager for fields over a mesh.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, mesh):
    """
    Constructor.
    """
    import pylith.topology.topology as bindings
    self.cppHandle = bindings.FieldsManager(mesh.cppHandle)
    return


  def addReal(self, label):
    """
    Create real field over mesh.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.addReal(label)


  def getReal(self, label):
    """
    Get real field over mesh.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.getReal(label)


  def delReal(self, label):
    """
    Delete real field over mesh.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.delReal(label)


  def setFiberDimension(self, label, fiberDim, points="vertices"):
    """
    Set fiber dimension for field

    points = { 'vertices', 'cells' }
    """
    assert(None != self.cppHandle)
    return self.cppHandle.setFiberDimension(label, fiberDim, points)


  def allocate(self, label):
    """
    Allocate field.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.allocate(label)


  def copyLayout(self, label):
    """
    Copy layout of field to all fields in manager.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.copyLayout(label)


  def copyLayoutFromSrc(self, field):
    """
    Copy layout of field to all fields in manager..
    """
    assert(None != self.cppHandle)
    return self.cppHandle.copyLayoutFromSrc(field)


  def solutionField(self, name):
    """
    Set name of field corresponding to solution.
    """
    assert(None != self.cppHandle)
    self.cppHandle.solutionField(name)
    return


  def getSolution(self):
    """
    Get field corresponding to solution.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.getSolution()


  def createHistory(self, labels):
    """
    Create history manager for a subset of the managed fields.
    """
    assert(None != self.cppHandle)
    self.cppHandle.createHistory(labels)
    return


  def shiftHistory(self):
    """
    Shift fields in history.
    """
    assert(None != self.cppHandle)
    self.cppHandle.shiftHistory()
    return


  def getHistoryItem(self, index):
    """
    Get field in history by position.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.getHistoryItem(index)


# End of file 
