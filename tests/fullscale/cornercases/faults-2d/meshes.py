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
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=16, ncorners=3, nvertices=14+2),

        # Materials
        "mat_xneg": MeshEntity(ncells=8, ncorners=3, nvertices=8),
        "mat_xpos": MeshEntity(ncells=8, ncorners=3, nvertices=8),

        # Faults
        "fault": MeshEntity(ncells=1, ncorners=2, nvertices=2),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_xpos": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_yneg": MeshEntity(ncells=4, ncorners=2, nvertices=5+1),
        "bc_ypos": MeshEntity(ncells=4, ncorners=2, nvertices=5+1),
    }


class Quad(object):
    """Mesh information for quad mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=4, ncorners=4, nvertices=10+2),

        # Materials
        "mat_xneg": MeshEntity(ncells=2, ncorners=4, nvertices=6),
        "mat_xpos": MeshEntity(ncells=2, ncorners=4, nvertices=6),

        # Faults
        "fault": MeshEntity(ncells=1, ncorners=2, nvertices=2),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_xpos": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_yneg": MeshEntity(ncells=4, ncorners=2, nvertices=5+1),
        "bc_ypos": MeshEntity(ncells=4, ncorners=2, nvertices=5+1),
    }


# End of file
