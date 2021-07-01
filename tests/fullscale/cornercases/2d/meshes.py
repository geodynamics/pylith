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
#
# @file tests/fullscale/cornercases/2d/meshes.py
#
# @brief Mesh information for test cases.


from pylith.testing.FullTestApp import MeshEntity


class Tri(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=2, ncorners=3, nvertices=4),

        # Materials
        "elastic": MeshEntity(ncells=2, ncorners=3, nvertices=4),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_xpos": MeshEntity(ncells=1, ncorners=2, nvertices=2),
    }


class Quad(object):
    """Mesh information for quad mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=1, ncorners=4, nvertices=4),

        # Materials
        "elastic": MeshEntity(ncells=1, ncorners=4, nvertices=4),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_xpos": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_yneg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_ypos": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "bc_domain": MeshEntity(ncells=4, ncorners=2, nvertices=4),
    }


# End of file
