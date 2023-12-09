# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.testing.FullTestApp import MeshEntity


class Tri(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=8, ncorners=3, nvertices=9+3),

        # Materials
        "elastic": MeshEntity(ncells=8, ncorners=3, nvertices=9+3),

        # Faults
        "fault": MeshEntity(ncells=2, ncorners=2, nvertices=3),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=2, ncorners=2, nvertices=3),
        "bc_xpos": MeshEntity(ncells=2, ncorners=2, nvertices=3),
        "bc_yneg": MeshEntity(ncells=2, ncorners=2, nvertices=3+1),
        "bc_ypos": MeshEntity(ncells=2, ncorners=2, nvertices=3+1),
    }


class Quad(object):
    """Mesh information for quad mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=4, ncorners=4, nvertices=9+3),

        # Materials
        "elastic": MeshEntity(ncells=4, ncorners=4, nvertices=9+3),

        # Faults
        "fault": MeshEntity(ncells=2, ncorners=2, nvertices=3),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=2, ncorners=2, nvertices=3),
        "bc_xpos": MeshEntity(ncells=2, ncorners=2, nvertices=3),
        "bc_yneg": MeshEntity(ncells=2, ncorners=2, nvertices=3+1),
        "bc_ypos": MeshEntity(ncells=2, ncorners=2, nvertices=3+1),
    }


# End of file
