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

##########  Somehow need to get filenames in here.
  def mesh(self):
    """Get a finite element mesh."""

    from lithomop3d.Mesh import Mesh
    mesh = Mesh()

    # Get nodes
    self._getNodes(mesh.nodes)

    # Get elements and define element families and materials
    self._getElements(mesh.elements)


    return mesh
    
  def __init__(self, name="meshimporter"):
    """Constructor."""
    Mesher.__init__(self, name)
    return
        

# version
__id__ = "$Id$"

# End of file 
