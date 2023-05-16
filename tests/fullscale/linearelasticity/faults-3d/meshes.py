#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

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
