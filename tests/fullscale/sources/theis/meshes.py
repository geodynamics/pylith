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
# @file tests/fullscale/poroelasticty/theis/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity


class Tri(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=200, ncorners=3, nvertices=202),

        # Materials
        "poroelastic": MeshEntity(ncells=200, ncorners=3, nvertices=202),

        # Boundaries
        "x_neg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "x_pos": MeshEntity(ncells=100, ncorners=2, nvertices=101),
        "y_neg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "y_pos": MeshEntity(ncells=100, ncorners=2, nvertices=101),
    }


class Quad(object):
    """Mesh information for quad mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=100, ncorners=4, nvertices=202),

        # Materials
        "poroelastic": MeshEntity(ncells=100, ncorners=4, nvertices=202),

        # Boundaries
        "x_neg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "x_pos": MeshEntity(ncells=100, ncorners=2, nvertices=101),
        "y_neg": MeshEntity(ncells=1, ncorners=2, nvertices=2),
        "y_pos": MeshEntity(ncells=100, ncorners=2, nvertices=101),
    }


# End of file
