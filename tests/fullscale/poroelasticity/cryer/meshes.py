# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/fullscale/poroelasticty/cryer/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity


class Tet(object):
    """
    Mesh information for tet mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=2735, ncorners=4, nvertices=718),

        # Materials
        "poroelastic": MeshEntity(ncells=2735, ncorners=4, nvertices=718),

        # Boundaries
        "bc_shell_traction": MeshEntity(ncells=404, ncorners=3, nvertices=227),
        "bc_shell_pressure": MeshEntity(ncells=404, ncorners=3, nvertices=227),
        "bc_xneg": MeshEntity(ncells=198, ncorners=3, nvertices=118),
        "bc_yneg": MeshEntity(ncells=196, ncorners=3, nvertices=117),
        "bc_zneg": MeshEntity(ncells=196, ncorners=3, nvertices=117),
    }


class Hex(object):
    """
    Mesh information for hex mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=896, ncorners=8, nvertices=1163),

        # Materials
        "poroelastic": MeshEntity(ncells=896, ncorners=8, nvertices=1163),

        # Boundaries
        "bc_shell_traction": MeshEntity(ncells=192, ncorners=4, nvertices=217),
        "bc_shell_pressure": MeshEntity(ncells=192, ncorners=4, nvertices=217),
        "x_neg": MeshEntity(ncells=96, ncorners=4, nvertices=115),
        "y_neg": MeshEntity(ncells=96, ncorners=4, nvertices=115),
        "z_neg": MeshEntity(ncells=96, ncorners=4, nvertices=115),
    }


# End of file
