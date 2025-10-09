# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
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
        "domain": MeshEntity(ncells=488, ncorners=3, nvertices=290),

        # Materials
        "poroelastic": MeshEntity(ncells=488, ncorners=3, nvertices=290),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=5, ncorners=2, nvertices=6),
        "bc_xpos": MeshEntity(ncells=5, ncorners=2, nvertices=6),
        "bc_yneg": MeshEntity(ncells=40, ncorners=2, nvertices=41),
        "bc_ypos": MeshEntity(ncells=40, ncorners=2, nvertices=41),
    }


class Quad(object):
    """Mesh information for quad mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=200, ncorners=4, nvertices=246),

        # Materials
        "poroelastic": MeshEntity(ncells=200, ncorners=4, nvertices=246),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=5, ncorners=2, nvertices=6),
        "bc_xpos": MeshEntity(ncells=5, ncorners=2, nvertices=6),
        "bc_yneg": MeshEntity(ncells=40, ncorners=2, nvertices=41),
        "bc_ypos": MeshEntity(ncells=40, ncorners=2, nvertices=41),
    }


# End of file
