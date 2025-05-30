# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/fullscale/linearelasticity/nofaults/2d/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity

class Tet(object):
    """Mesh information for tet mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=374, ncorners=4, nvertices=112),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "upper_crust": MeshEntity(ncells=140, ncorners=4, nvertices=57),
        "lower_crust": MeshEntity(ncells=234, ncorners=4, nvertices=80),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=27, ncorners=3, nvertices=19),
        "bc_xpos": MeshEntity(ncells=27, ncorners=3, nvertices=19),
        "bc_yneg": MeshEntity(ncells=27, ncorners=3, nvertices=19),
        "bc_ypos": MeshEntity(ncells=27, ncorners=3, nvertices=19),
        "bc_zneg": MeshEntity(ncells=32, ncorners=3, nvertices=25),
        "groundsurf": MeshEntity(ncells=32, ncorners=3, nvertices=25),
    }


class Hex(object):
    """Mesh information for hex mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=64, ncorners=8, nvertices=100),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "upper_crust": MeshEntity(ncells=16, ncorners=8, nvertices=50),
        "lower_crust": MeshEntity(ncells=32, ncorners=8, nvertices=75),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=12, ncorners=4, nvertices=20),
        "bc_xpos": MeshEntity(ncells=12, ncorners=4, nvertices=20),
        "bc_yneg": MeshEntity(ncells=12, ncorners=4, nvertices=20),
        "bc_ypos": MeshEntity(ncells=12, ncorners=4, nvertices=20),
        "bc_zneg": MeshEntity(ncells=16, ncorners=4, nvertices=25),
        "groundsurf": MeshEntity(ncells=16, ncorners=4, nvertices=25),
    }


# End of file
