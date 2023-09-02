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
    """Mesh information for tet mesh using Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=1450, ncorners=4, nvertices=425+9),

        # Materials
        "mat_elastic": MeshEntity(ncells=1450, ncorners=4, nvertices=425+9),

        # Faults
        "fault": MeshEntity(ncells=22, ncorners=3, nvertices=18),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=84, ncorners=3, nvertices=55),
        "bc_xpos": MeshEntity(ncells=84, ncorners=3, nvertices=55),
        "bc_yneg": MeshEntity(ncells=84, ncorners=3, nvertices=55),
        "bc_ypos": MeshEntity(ncells=84, ncorners=3, nvertices=55),
        "bc_zneg": MeshEntity(ncells=162, ncorners=3, nvertices=98),
        "bc_zpos": MeshEntity(ncells=164, ncorners=3, nvertices=94),
    }


# End of file
