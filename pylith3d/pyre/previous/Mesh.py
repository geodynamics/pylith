#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
#
#  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
#  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
#  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Mesh class
class Mesh(object):
    """Template for a mesh object."""

    def __init__(self):
        """Constructor."""
        # This object will contain all of the objects below.
        return
        
# Nodes class
class Nodes(object):
    """Template for an object holding nodal coordinate info."""

    def __init__(self):
        """Constructor."""
        # self.numberSpaceDimensions = 0
        self.numberNodes = 0
        self.coordArray = [0.0]
        self.pointerToCoordArray = None
        return

# Elements class
class Elements(object):
    """Template for an object holding element info."""

    def __init__(self):
        """Constructor."""
        self.numberElements = 0
        # This object can contain an ElementNodes class and, optionally,
        # a Connectivities class, although I don't need that.
        return

# ElementNodes class
class ElementNodes(object):
    """Template for an ElementNodes object."""

    def __init__(self):
        self.elementType = 0
        self.elementNodeList = [0]
        self.pointerToElementNodeList = None
        return

# Connectivities class
class Connectivities(object):
    """Template for a Connectivities object."""

    def __init__(self):
        self.numberElementNeighbors = 0
        self.elementNeighborsList = [0]
        self.pointerToElementNeighborsList = None
        return

# NodalBoundaries class:
class NodalBoundaries(object):
    """Template for a NodalBoundaries object."""

    def __init__(self):
        self.numberNodalBoundaries = 0
        # This object will contain the NodalBoundary object (for a single boundary)
        # described below.
        return

# NodalBoundary class:
class NodalBoundary(object):
    """Template for a NodalBoundary class."""

    def __init__(self):
        self.numberBoundaryNodes = 0
        self.boundaryNodesList = [0]
        self.pointerToBoundaryNodesList = None

# ElementBoundaries class:
class ElementBoundaries(object):
    """Template for an ElementBoundaries class."""

    def __init__(self):
        self.numberElementBoundaries = 0
        # This object will contain the ElementBoundary object (for a single boundary)
        # described below.
        return

# ElementBoundary class
class ElementBoundary(object):
    """Template for an ElementBoundary object."""

    def __init__(self):
        self.numberBoundaryFacets = 0
        # This object will contain the BoundaryElementNodes object (for a single
        # facet) described below.
        return

# BoundaryElementNodes class
class BoundaryElementNodes(object):
    """Template for a BoundaryElementNodes object."""

    def __init__(self):
        self.facetType = 0
        self.facetNodeList = [0]
        self.pointerToFacetNodeList = None
        return
    
# version
__id__ = "$Id: Mesh.py,v 1.1 2004/10/11 15:50:26 willic3 Exp $"

# Generated automatically by PythonMill on Tue Oct  5 14:06:12 2004

# End of file 
