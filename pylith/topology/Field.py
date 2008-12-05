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

## @file pylith/topology/Field.py
##
## @brief Python object for managing a vector field over vertices or
## cells of a finite-element mesh.
##
## Factory: vector_field

# Field class
class Field(object):
  """
  Python object for managing a vector field over vertices or cells of
  a finite-element mesh.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, mesh):
    """
    Constructor.
    """
    self.cppHandle = self._createCppHandle(mesh)
    return


  def setName(self, value):
    """
    Set name of field.
    """
    assert(None != self.cppHandle)
    self.cppHandle.name = value
    return


  def getName(self):
    """
    Get name of field.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.name


  def setVectorFieldType(self, value):
    """
    Set vector field type.
    """
    assert(None != self.cppHandle)
    self.cppHandle.vectorFieldType = value
    return


  def getVectorFieldType(self):
    """
    Get vector field type.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.VectorFieldType
    return


  def setSpaceDim(self, value):
    """
    Set spaceDim of field.
    """
    assert(None != self.cppHandle)
    self.cppHandle.spaceDim = value
    return


  def getSpaceDim(self, spaceDim):
    """
    Get spaceDim of field.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.spaceDim


  def setScale(self, value):
    """
    Set scale of field.
    """
    assert(None != self.cppHandle)
    self.cppHandle.scale = value
    return


  def getScale(self):
    """
    Get scale of field.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.scale


  def setAddDimensionOkay(self, value):
    """
    Set addDimensionOkay.
    """
    assert(None != self.cppHandle)
    self.cppHandle.addDimensionOkay = value
    return


  def getAddDimensionOkay(self):
    """
    Get addDimensionOkay.
    """
    assert(None != self.cppHandle)
    return self.cppHandle.addDimensionOkay


  def copyLayout(self, field):
    """
    Copy layout of field.
    """
    assert(None != self.cppHandle)
    self.cppHandle.copyLayout(field.cppHandle)
    return


  def clear(self):
    """
    Clear variables associated with section.
    """
    assert(None != self.cppHandle)
    self.cppHandle.clear()
    return


  def zero(self):
    """
    Zero section values.
    """
    assert(None != self.cppHandle)
    self.cppHandle.zero()
    return
    
    
  def complete(self):
    """
    Complete section by assembling over processors.
    """
    assert(None != self.cppHandle)
    self.cppHandle.complete()
    return
    
    
  def copy(self, field):
    """
    Copy field values and metadata.
    """
    assert(None != self.cppHandle)
    self.cppHandle.copy(field.cppHandle)
    return
    
    
  def add(self, field):
    """
    Add two fields, storing result in field.
    """
    assert(None != self.cppHandle)
    self.cppHandle.add(field.cppHandle)
    return
    

  def dimensionalize(self):
    """
    Add dimensions to field.
    """
    okay = self.cppHandle.addDimensionOkay
    if ~okay:
      name = self.cppHandle.name
      raise RuntimeError("Field '%s' is protected. Cannot dimensionalize." % \
                         name)
    assert(None != self.cppHandle)
    self.cppHandle.dimensionalize()
    

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.topology.topology as bindings
      self.cppHandle = bindings.Field()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def vector_field():
  """
  Factory associated with Mesh.
  """
  return Field()


# End of file
