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
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

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
