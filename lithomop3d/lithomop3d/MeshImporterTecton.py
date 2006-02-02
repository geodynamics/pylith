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
  """Class for importing a mesh using Tecton-style (original LithoMop) input files."""

  class Inventory(MeshImporter.Inventory):

    import pyre.inventory

    coordFile = pyre.inventory.str("coordFile",default="")
    elemFile = pyre.inventory.str("elemFile",default="")
    f77FileInput = pyre.inventory.int("f77FileInput", default=10)


  def generate(self, fileRoot):
    """Get a finite element mesh."""

    from lithomop3d.Mesh import Mesh
    mesh = Mesh()

    # Get nodes
    mesh.nodes = self._getNodes(fileRoot)

    # Get elements and define element families and materials
    mesh.elements, mesh.elemFamily = self._getElems(fileRoot)

    return mesh

  def _getNodes(self, fileRoot):
    """Gets dimensions and nodal coordinates from Tecton-style input file."""

    import pyre.units
    import lithomop3d

    # Get input filename
    if self.inventory.coordFile == "":
      coordInputFile = fileRoot + ".coord"
    else:
      coordInputFile = self.inventory.coordFile

    f77FileInput = self.inventory.f77FileInput
    
    numNodes, numDims, coordUnits = lithomop3d.scan_coords(
      f77FileInput,
      coordInputFile)

    coordScaleString = \
                     uparser.parse(string.strip(coordUnits))
    coordScaleFactor = coordScaleString.value

    ptrCoords = lithomop3d.allocateDouble(numDims*numNodes)

    nodes = { 'numNodes': numNodes,
              'dim': numDims,
              'ptrCoords': ptrCoords }

    lithomop3d.read_coords(nodes,
                           coordScaleFactor,
                           f77FileInput,
                           coordInputFile)

    return nodes
      

  def _getElems(self, fileRoot):
    """Gets dimensions and element nodes from Tecton-style input file, and
    defines element families based on element group number."""

    import lithomop3d

    # Get input filename
    if self.inventory.elemFile == "":
      elemInputFile = fileRoot + ".connect"
    else:
      elemInputFile = self.inventory.elemFile

    f77FileInput = self.inventory.f77FileInput
    

    ptrNumNodesPerElemType = lithomop3d.intListToArray(mesh.numNodesPerElemType)
    numElemTypes = len(numNodesPerElemType)
    maxNumElemFamilies = 1024
    ptrTmpElemFamilySizes = lithomop3d.allocateInt(maxNumElemFamilies)

    numElems, numNodesPerElem, elemType, numElemFamilies = lithomop3d.scan_connect(
      ptrNumNodesPerElemType,
      numElemTypes,
      ptrTmpElemFamilySizes,
      maxNumElemFamilies,
      f77FileInput,
      elemInputFile)

    ptrConns = lithomop3d.allocateInt(numElems*numNodesPerElem)
    ptrInitOrder = lithomop3d.allocateInt(numElems)

    elements = {'numElems': numElems,
                'numNodesPerElem': numNodesPerElem,
                'elemType': elemType,
                'ptrConns': ptrConns,
                'ptrInitOrder': ptrInitOrder}

    ptrElemFamilySizes = lithomop3d.allocateInt(numElemFamilies)

    elemFamily = {'numElemFamilies': numElemFamilies,
                  'ptrElemFamilySizes': ptrElemFamilySizes}
      
    lithomop3d.read_connect(
      elements,
      elemFamily,
      ptrTmpElemFamilySizes,
      maxNumElemFamilies,
      f77FileInput,
      elemInputFile)

    ptrTmpElemFamilySizes = None

    return elements, elemFamily
      
    
  def __init__(self, name="meshimporter"):
    """Constructor."""

    Mesher.__init__(self, name)

    return
        

# version
__id__ = "$Id$"

# End of file 
