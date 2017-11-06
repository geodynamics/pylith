#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pyre/meshio/VertexFilterDecimateP1.py
##
# @brief Python class for decimating vertex field to P1.
##
# Factory: output_vertex_filter

from VertexFilter import VertexFilter
from meshio import VertexFilterDecimateP1 as ModuleVertexFilterDecimateP1

# VertexFilterDecimateP1 class


class VertexFilterDecimateP1(VertexFilter, ModuleVertexFilterDecimateP1):
    """
    Python class for decimating vertex field to P1.

    Factory: output_vertex_filter
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="vertexfilterdecimatep1"):
        """
        Constructor.
        """
        VertexFilter.__init__(self, name)
        ModuleVertexFilterDecimateP1.__init__(self)
        self.filter = True
        return


# FACTORIES ////////////////////////////////////////////////////////////

def output_vertex_filter():
    """
    Factory associated with VertexFilter.
    """
    return VertexFilterDecimateP1()


# End of file
