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
# @file tests/fullscale/poroelasticty/cryer/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity


class Tet(object):
    """
    Mesh information for tet mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=5398, ncorners=4, nvertices=1150),

        # Materials
        "poroelastic": MeshEntity(ncells=5398, ncorners=4, nvertices=1150),

        # Boundaries
        "surface_traction": MeshEntity(ncells=374, ncorners=3, nvertices=212),
        "surface_pressure": MeshEntity(ncells=374, ncorners=3, nvertices=212),
        "x_neg": MeshEntity(ncells=184, ncorners=3, nvertices=111),
        "y_neg": MeshEntity(ncells=184, ncorners=3, nvertices=111),
        "z_neg": MeshEntity(ncells=184, ncorners=3, nvertices=111),
    }


class Hex(object):
    """
    Mesh information for hex mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=896, ncorners=8, nvertices=1163),

        # Materials
        "poroelastic": MeshEntity(ncells=896, ncorners=8, nvertices=1163),

        # Boundaries
        "surface_traction": MeshEntity(ncells=192, ncorners=4, nvertices=217),
        "surface_pressure": MeshEntity(ncells=192, ncorners=4, nvertices=217),
        "x_neg": MeshEntity(ncells=96, ncorners=4, nvertices=115),
        "y_neg": MeshEntity(ncells=96, ncorners=4, nvertices=115),
        "z_neg": MeshEntity(ncells=96, ncorners=4, nvertices=115),
    }


# End of file
