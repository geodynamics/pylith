#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#             Copyright (C) 2006 Rensselaer Polytechnic Institute
#
# 
# 	Permission is hereby granted, free of charge, to any person
# 	obtaining a copy of this software and associated documentation
# 	files (the "Software"), to deal in the Software without
# 	restriction, including without limitation the rights to use,
# 	copy, modify, merge, publish, distribute, sublicense, and/or
# 	sell copies of the Software, and to permit persons to whom the
# 	Software is furnished to do so, subject to the following
# 	conditions:
# 
# 	The above copyright notice and this permission notice shall be
# 	included in all copies or substantial portions of the Software.
# 
# 	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# 	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# 	OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# 	NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# 	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# 	WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# 	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# 	OTHER DEALINGS IN THE SOFTWARE.
#         
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Initial attempt at a module that gets a Mesh object from Tecton-style input files.

class MeshImporterTecton(MeshImporter):
  """Class for importing a mesh using Tecton-style (original PyLith) input files."""

  class Inventory(MeshImporter.Inventory):

    import pyre.inventory

    coordFile = pyre.inventory.str("coordFile",default="")
    elemFile = pyre.inventory.str("elemFile",default="")
    f77FileInput = pyre.inventory.int("f77FileInput", default=10)
    maxNumElemFamilies = pyre.inventory.int("maxNumElemFamilies", default=1024)


  def generate(self, fileRoot):
    """Get a finite element mesh."""

    from pylith3d.Mesh import Mesh
    mesh = Mesh()

    # Get nodes
    mesh.nodes = self._getNodes(fileRoot)

    # Get elements and define element families and materials
    mesh.elements, mesh.elemFamilies = self._getElems(fileRoot)

    return mesh

  def _getNodes(self, fileRoot):
    """Gets dimensions and nodal coordinates from Tecton-style input file."""

    import pyre.units
    import pylith3d as pl3d

    # Get input filename
    if self.inventory.coordFile == "":
      coordInputFile = fileRoot + ".coord"
    else:
      coordInputFile = self.inventory.coordFile

    f77FileInput = self.inventory.f77FileInput
    
    numNodes, numDims, coordUnits = pl3d.scan_coords(
      f77FileInput,
      coordInputFile)

    coordScaleString = \
                     uparser.parse(string.strip(coordUnits))
    coordScaleFactor = coordScaleString.value

    ptrCoords = pl3d.allocateDouble(numDims*numNodes)

    nodes = { 'numNodes': numNodes,
              'dim': numDims,
              'ptrCoords': ptrCoords }

    pl3d.read_coords(nodes['ptrCoords'],
                     coordScaleFactor,
                     nodes['numNodes'],
                     nodes['dim'],
                     f77FileInput,
                     coordInputFile)

    return nodes
      

  def _getElems(self, fileRoot):
    """Gets dimensions and element nodes from Tecton-style input file, and
    defines element families based on element group number."""

    import pylith3d as pl3d

    # Get input filename
    if self.inventory.elemFile == "":
      elemInputFile = fileRoot + ".connect"
    else:
      elemInputFile = self.inventory.elemFile

    f77FileInput = self.inventory.f77FileInput
    maxNumElemFamilies = self.inventory.maxNumElemFamilies
    
    numElemTypes = len(mesh.numNodesPerElemType)
    # I am changing the way this is done for now, under the assumption that we
    # will be using f2py to generate bindings.  This will allow me to 'see'
    # specified arrays as lists in python.
    # ptrTmpElemFamilySizes = pylith3d.allocateInt(maxNumElemFamilies)
    tmpElemFamilySizes = [0]*maxNumElemFamilies
    tmpElemFamilyTypes = [0]*maxNumElemFamilies
    tmpElemFamilyMatIds = [0]*maxNumElemFamilies

    # I need to check that I've included everything here, and I will also need
    # to change the routine arguments and bindings.
    numElems,
    numElemFamilies,
    tmpElemFamilySizes,
    tmpElemFamilyTypes,
    tmpElemFamilyMatIds = pl3d.scan_connect(
      mesh.numNodesPerElemType,
      numElemTypes,
      maxNumElemFamilies,
      f77FileInput,
      elemInputFile)

    # Compact temporary element family arrays to get actual arrays, and
    # determine size of element node array.  Store information in a list of
    # dictionaries.
    elemFamilies = []
    elemNodeArraySize = 0
    count = 0

    for i in range(maxNumElemFamilies):
      if tmpElemFamilySizes[i] != 0 and tmpElemFamilyTypes[i] != 0 and tmpElemFamilyMatIds[i] != 0:
        elemFamilyNodesSize = tmpElemFamilySizes[i]*mesh.numNodesPerElemType[tmpElemFamilyTypes[i]]
        elemNodeArraySize += elemFamilyNodesSize
        ptrElemFamilyNodes = pl3d.allocateInt(elemFamilyNodesSize)
        elemFamilies.append({'elemFamilySize': tmpElemFamilySizes[i],
                             'elemFamilyType': tmpElemFamilyTypes[i],
                             'elemFamilyMatId': tmpElemFamilyMatIds[i],
                             'elemFamilyName': "family" + str(count),
                             'ptrElemFamilyNodes': ptrElemFamilyNodes})
        count=count + 1

    tmpElemFamilySizes = None
    tmpElemFamilyTypes = None
    tmpElemFamilyMatIds = None

    elemFamilyError = "ElemFamily Error"
    if len(elemFamilies) != numElemFamilies:
      raise elemFamilyError, "Number of element families does not match input!"

    # Read the element node array.  The way things are set up now, I need to provide
    # arrays to contain the material ID and element type for each element
    # until the elements are sorted.
    ptrElemNodeArray = pl3d.allocateInt(elemNodeArraySize)
    ptrElemTypes = pl3d.allocateInt(numElems)
    ptrElemIds = pl3d.allocateInt(numElems)

    pylith3d.read_connect(
      ptrElemNodeArray,
      ptrElemTypes,
      ptrElemIds,
      mesh.numNodesPerElemType,
      elemNodeArraySize,
      numElems,
      numElemTypes,
      f77FileInput,
      elemInputFile)

    # Create elements dictionary.  This info is no longer needed after element
    # families have been formed.
    elements = {'numElems': numElems,
                'elemNodeArraySize': elemNodeArraySize,
                'ptrElemNodeArray': ptrElemNodeArray,
                'ptrElemTypes': ptrElemTypes,
                'ptrElemIds': ptrElemIds}

    return elements, elemFamilies
      
    
  def __init__(self, name="meshimporter"):
    """Constructor."""

    Mesher.__init__(self, name)

    return
        

# version
__id__ = "$Id$"

# End of file 
