# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pylith.testing.FullTestApp import MeshEntity


class Tri(object):
    """Mesh information for tri mesh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=206, ncorners=3, nvertices=156),
        # Materials
        "poroelastic": MeshEntity(ncells=206, ncorners=3, nvertices=156),
        # Boundaries
        "bc_xneg": MeshEntity(ncells=50, ncorners=2, nvertices=51),
        "bc_xpos": MeshEntity(ncells=50, ncorners=2, nvertices=51),
        "bc_yneg": MeshEntity(ncells=2, ncorners=2, nvertices=3),
        "bc_ypos_pressure": MeshEntity(ncells=2, ncorners=2, nvertices=3),
        "bc_ypos_traction": MeshEntity(ncells=2, ncorners=2, nvertices=3),
    }


class Quad(object):
    """Mesh information for quad mesh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=100, ncorners=4, nvertices=153),
        # Materials
        "poroelastic": MeshEntity(ncells=100, ncorners=4, nvertices=153),
        # Boundaries
        "bc_xneg": MeshEntity(ncells=50, ncorners=2, nvertices=51),
        "bc_xpos": MeshEntity(ncells=50, ncorners=2, nvertices=51),
        "bc_yneg": MeshEntity(ncells=2, ncorners=2, nvertices=3),
        "bc_ypos_pressure": MeshEntity(ncells=2, ncorners=2, nvertices=3),
        "bc_ypos_traction": MeshEntity(ncells=2, ncorners=2, nvertices=3),
    }


# End of file
