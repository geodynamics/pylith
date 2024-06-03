# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/fullscale/cornercases/3d/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity


class Tet(object):
    """Mesh information for tet mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=5, ncorners=4, nvertices=8),

        # Materials
        "elastic": MeshEntity(ncells=5, ncorners=4, nvertices=8),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=2, ncorners=3, nvertices=4),
        "bc_xpos": MeshEntity(ncells=2, ncorners=3, nvertices=4),
    }


class Hex(object):
    """Mesh information for hex mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=1, ncorners=8, nvertices=8),

        # Materials
        "elastic": MeshEntity(ncells=1, ncorners=8, nvertices=8),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=1, ncorners=4, nvertices=4),
        "bc_xpos": MeshEntity(ncells=1, ncorners=4, nvertices=4),
        "bc_yneg": MeshEntity(ncells=1, ncorners=4, nvertices=4),
        "bc_ypos": MeshEntity(ncells=1, ncorners=4, nvertices=4),
        "bc_zneg": MeshEntity(ncells=1, ncorners=4, nvertices=4),
        "bc_zpos": MeshEntity(ncells=1, ncorners=4, nvertices=4),
        "bc_domain": MeshEntity(ncells=6, ncorners=4, nvertices=8),
    }


# End of file
