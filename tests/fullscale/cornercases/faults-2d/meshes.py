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
