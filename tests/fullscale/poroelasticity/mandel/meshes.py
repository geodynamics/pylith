# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/fullscale/poroelasticty/mandel/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity


class Tri(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=320, ncorners=3, nvertices=205),

        # Materials
        "poroelastic": MeshEntity(ncells=320, ncorners=3, nvertices=205),

        # Boundaries
        "x_neg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "x_pos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "y_neg": MeshEntity(ncells=40, ncorners=2, nvertices=41),
        "y_pos": MeshEntity(ncells=40, ncorners=2, nvertices=41),
    }


class Quad(object):
    """Mesh information for quad mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=160, ncorners=4, nvertices=205),

        # Materials
        "poroelastic": MeshEntity(ncells=160, ncorners=4, nvertices=205),

        # Boundaries
        "x_neg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "x_pos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "y_neg": MeshEntity(ncells=40, ncorners=2, nvertices=41),
        "y_pos": MeshEntity(ncells=40, ncorners=2, nvertices=41),
    }


# End of file
