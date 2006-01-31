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

    coordFile = pyre.inventory.str("coordFile")
    elemFile = pyre.inventory.str("elemFile")
    f77FileInput = pyre.inventory.int("f77FileInput", default=10)


  def generate(self):
    """Get a finite element mesh."""

    from lithomop3d.Mesh import Mesh
    mesh = Mesh()

    # Get nodes
    mesh.nodes = self._getNodes()

    # Get elements and define element families and materials
    mesh.elements, mesh.elemFamily = self._getElems(mesh.numNodesPerElemType)


    return mesh

  def _getNodes(self):
    """Gets dimensions and nodal coordinates from Tecton-style input file."""

    import pyre.units
    import lithomop3d
    
    nodeInfo = lithomop3d.scan_coords(
      f77FileInput,
      numDims,
      coordInputFile)

    numNodes = nodeInfo[0]
    numDims = nodeInfo[1]
    coordUnits = nodeInfo[2]
    
    coordScaleString = \
                     uparser.parse(string.strip(coordUnits))
    coordScaleFactor = coordScaleString.value

    ptrCoords = lithomop3d.allocateDouble(numDims*numNodes)

    nodes['numNodes'] = numNodes
    nodes['dim'] = numDims
    nodes['ptrCoords'] = ptrCoords

    lithomop3d.read_coords(nodes,
                           coordScaleFactor,
                           self.inventory.f77FileInput,
                           self.inventory.coordInputFile)

    return nodes
      
  def _getElems(self, numNodesPerElemType):
    """Gets dimensions and element nodes from Tecton-style input file, and
    defines element families based on material model."""

    import lithomop3d

    ptrNumNodesPerElemType = lithomop3d.intListToArray(numNodesPerElemType)
    numElemTypes = len(numNodesPerElemType)

    elemInfo = lithomop3d.scan_connect(
      ptrNumNodesPerElemType,
      
      

    
  def __init__(self, name="meshimporter"):
    """Constructor."""

    Mesher.__init__(self, name)

    return
        

# version
__id__ = "$Id$"

# End of file 
