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
        "domain": MeshEntity(ncells=44, ncorners=3, nvertices=31 + 5),
        # Materials
        "poroelastic": MeshEntity(ncells=44, ncorners=3, nvertices=31 + 5),
        # Faults
        "fault": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        # Boundaries
        "bc_disp_xneg": MeshEntity(ncells=4, ncorners=2, nvertices=5 + 1),
        "bc_disp_xpos": MeshEntity(ncells=4, ncorners=2, nvertices=5 + 1),
        "bc_disp_yneg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_disp_ypos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_press_xneg": MeshEntity(ncells=4, ncorners=2, nvertices=5 + 1),
        "bc_press_xpos": MeshEntity(ncells=4, ncorners=2, nvertices=5 + 1),
    }


class QuadGmsh(object):
    """Mesh information for quad mesh using Gmsh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=16, ncorners=4, nvertices=25 + 5),
        # Materials
        "poroelastic": MeshEntity(ncells=16, ncorners=3, nvertices=25 + 5),
        # Faults
        "fault": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        # Boundaries
        "bc_disp_xneg": MeshEntity(ncells=4, ncorners=2, nvertices=5 + 1),
        "bc_disp_xpos": MeshEntity(ncells=4, ncorners=2, nvertices=5 + 1),
        "bc_disp_yneg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_disp_ypos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_press_xneg": MeshEntity(ncells=4, ncorners=2, nvertices=5 + 1),
        "bc_press_xpos": MeshEntity(ncells=4, ncorners=2, nvertices=5 + 1),
    }


# End of file
