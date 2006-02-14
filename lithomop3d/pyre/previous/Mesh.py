#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2004  All Rights Reserved
#
#  Copyright 2004 Rensselaer Polytechnic Institute.
#  All worldwide rights reserved.  A license to use, copy, modify and
#  distribute this software for non-commercial research purposes only
#  is hereby granted, provided that this copyright notice and
#  accompanying disclaimer is not modified or removed from the software.
#
#  DISCLAIMER:  The software is distributed "AS IS" without any express
#  or implied warranty, including but not limited to, any implied
#  warranties of merchantability or fitness for a particular purpose
#  or any warranty of non-infringement of any current or pending patent
#  rights.  The authors of the software make no representations about
#  the suitability of this software for any particular purpose.  The
#  entire risk as to the quality and performance of the software is with
#  the user.  Should the software prove defective, the user assumes the
#  cost of all necessary servicing, repair or correction.  In
#  particular, neither Rensselaer Polytechnic Institute, nor the authors
#  of the software are liable for any indirect, special, consequential,
#  or incidental damages related to the software, to the maximum extent
#  the law permits.
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
