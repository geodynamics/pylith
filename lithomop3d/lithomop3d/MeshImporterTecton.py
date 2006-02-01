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

    coordFile = pyre.inventory.str("coordFile",default="None")
    elemFile = pyre.inventory.str("elemFile",default="None")
    f77FileInput = pyre.inventory.int("f77FileInput", default=10)


  def generate(self, fileRoot, ptrMatModInfo):
    """Get a finite element mesh."""

    from lithomop3d.Mesh import Mesh
    mesh = Mesh()

    # Get nodes
    mesh.nodes = self._getNodes(fileRoot)

    # Get elements and define element families and materials
    mesh.elements, mesh.elemFamily = self._getElems(fileRoot, ptrMatModInfo)


    return mesh

  def _getNodes(self, fileRoot):
    """Gets dimensions and nodal coordinates from Tecton-style input file."""

    import pyre.units
    import lithomop3d

    # Get input filename
    if self.inventory.coordFile == "None":
      coordInputFile = fileRoot + ".ascii"
    else:
      coordInputFile = self.inventory.coordFile
    
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
      
  def _getElems(self, fileRoot, ptrMatModInfo):
    """Gets dimensions and element nodes from Tecton-style input file, and
    defines element families based on material model."""

    import lithomop3d

    ptrNumNodesPerElemType = lithomop3d.intListToArray(mesh.numNodesPerElemType)
    numElemTypes = len(numNodesPerElemType)

    #  I need to figure out how to get most of the material stuff out of here.
    #  Current idea:
    #  1.  In the element nodes file, specify a materialGroup for each element.
    #  2.  In a separate file, identify a materialModel and spatialDatabase to be
    #      used for each group.
    #  3.  Ideally, the file identifying the model and database for each group would
    #      be read before the mesh info, so that I would already know how many groups
    #      there are.  If I don't do this, I won't know how large to allocate the array
    #      that keeps track of the number of elements in each group.
    #  For now, the number of groups could correspond to the number of element families,
    #  although I have to decide if this is better than having it correspond to the
    #  number of material models used (which will always be less than or equal to the
    #  number of groups).  In the future, the element families could be extended to also
    #  include element type.
    elemInfo = lithomop3d.scan_connect(
      ptrNumNodesPerElemType,
      
      

    
  def __init__(self, name="meshimporter"):
    """Constructor."""

    Mesher.__init__(self, name)

    return
        

# version
__id__ = "$Id$"

# End of file 
