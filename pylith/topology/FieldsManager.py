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
    return self.cppHandle.copyLayoutFromSrc(label)


# End of file 
