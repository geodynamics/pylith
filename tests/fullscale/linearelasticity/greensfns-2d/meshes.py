# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.testing.FullTestApp import MeshEntity


class TriGmsh(object):
    """Mesh information for tri mesh using Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=128, ncorners=3, nvertices=81+9),

        # Materials
        "mat_xneg": MeshEntity(ncells=64, ncorners=3, nvertices=45),
        "mat_xpos": MeshEntity(ncells=64, ncorners=3, nvertices=45),

        # Faults
        "fault": MeshEntity(ncells=8, ncorners=2, nvertices=9),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=10),
    }


class QuadGmsh(object):
    """Mesh information for quad mesh using Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=90, ncorners=4, nvertices=110+10),

        # Materials
        "mat_xneg": MeshEntity(ncells=45, ncorners=3, nvertices=60),
        "mat_xpos": MeshEntity(ncells=45, ncorners=3, nvertices=60),

        # Faults
        "fault": MeshEntity(ncells=9, ncorners=2, nvertices=10),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=9, ncorners=2, nvertices=10),
        "bc_xpos": MeshEntity(ncells=9, ncorners=2, nvertices=10),
        "bc_yneg": MeshEntity(ncells=11, ncorners=2, nvertices=12),
        "bc_ypos": MeshEntity(ncells=11, ncorners=2, nvertices=12),
    }


# End of file
