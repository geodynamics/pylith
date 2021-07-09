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
# @file tests/fullscale/poroelasticty/terzaghi/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity


class Tri(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=320, ncorners=3, nvertices=205),

        # Materials
        "poroelastic": MeshEntity(ncells=320, ncorners=3, nvertices=205),

        # Boundaries
        "x_neg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "x_pos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "y_neg": MeshEntity(ncells=40, ncorners=2, nvertices=41),
        "y_pos": MeshEntity(ncells=40, ncorners=2, nvertices=41),
    }


class Quad(object):
    """Mesh information for quad mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=160, ncorners=4, nvertices=205),

        # Materials
        "poroelastic": MeshEntity(ncells=160, ncorners=4, nvertices=205),

        # Boundaries
        "x_neg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "x_pos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "y_neg": MeshEntity(ncells=40, ncorners=2, nvertices=41),
        "y_pos": MeshEntity(ncells=40, ncorners=2, nvertices=41),
    }


# End of file
