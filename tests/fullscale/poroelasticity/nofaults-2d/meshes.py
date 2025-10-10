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


class TriGmsh(object):
    """Mesh information for tri mesh using Gmsh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=42, ncorners=3, nvertices=30),
        # Materials
        "poroelastic": MeshEntity(ncells=42, ncorners=3, nvertices=30),
        # Boundaries
        "bc_disp_xneg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_disp_xpos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_disp_yneg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_disp_ypos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_press_ypos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_press_xneg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
    }


class QuadGmsh(object):
    """Mesh information for quad mesh using Gmsh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=16, ncorners=4, nvertices=25),
        # Materials
        "poroelastic": MeshEntity(ncells=16, ncorners=3, nvertices=25),
        # Boundaries
        "bc_disp_xneg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_disp_xpos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_disp_yneg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_disp_ypos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_press_ypos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_press_xneg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
    }


# End of file
