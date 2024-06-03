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


class TetGmsh(object):
    """Mesh information for tri mesh using Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=571, ncorners=4, nvertices=182+38),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "mat_elastic": MeshEntity(ncells=571, ncorners=4, nvertices=182+38),

        # Faults
        "fault": MeshEntity(ncells=56, ncorners=3, nvertices=38),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=44, ncorners=3, nvertices=31),
        "bc_xpos": MeshEntity(ncells=44, ncorners=3, nvertices=31),
        "bc_yneg": MeshEntity(ncells=56, ncorners=3, nvertices=38+6),
        "bc_ypos": MeshEntity(ncells=56, ncorners=3, nvertices=38+6),
        "bc_zneg": MeshEntity(ncells=50, ncorners=3, nvertices=35+5),
        "bc_zpos": MeshEntity(ncells=50, ncorners=3, nvertices=35+5),
    }


class HexGmsh(object):
    """Mesh information for quad mesh using Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=150, ncorners=8, nvertices=252+36),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "mat_elastic": MeshEntity(ncells=150, ncorners=8, nvertices=252+36),

        # Faults
        "fault": MeshEntity(ncells=25, ncorners=4, nvertices=36),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=25, ncorners=4, nvertices=36),
        "bc_xpos": MeshEntity(ncells=25, ncorners=4, nvertices=36),
        "bc_yneg": MeshEntity(ncells=42, ncorners=4, nvertices=30+6),
        "bc_ypos": MeshEntity(ncells=42, ncorners=4, nvertices=30+6),
        "bc_zneg": MeshEntity(ncells=42, ncorners=4, nvertices=30+6),
        "bc_zpos": MeshEntity(ncells=42, ncorners=4, nvertices=30+6),
    }


# End of file
