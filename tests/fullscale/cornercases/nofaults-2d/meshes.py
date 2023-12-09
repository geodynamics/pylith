# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/fullscale/cornercases/2d/meshes.py
#
# @brief Mesh information for test cases.


from pylith.testing.FullTestApp import MeshEntity


class Tri(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=2, ncorners=3, nvertices=4),

        # Materials
        "elastic": MeshEntity(ncells=2, ncorners=3, nvertices=4),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_xpos": MeshEntity(ncells=1, ncorners=2, nvertices=2),
    }


class Quad(object):
    """Mesh information for quad mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=1, ncorners=4, nvertices=4),

        # Materials
        "elastic": MeshEntity(ncells=1, ncorners=4, nvertices=4),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_xpos": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_yneg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_ypos": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_domain": MeshEntity(ncells=4, ncorners=2, nvertices=4),
    }


# End of file
