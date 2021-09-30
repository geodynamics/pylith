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
# @file tests/fullscale/linearelasticity/cylinder-3d/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity

class Tet(object):
    """Mesh information for tet mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=10234, ncorners=4, nvertices=2361),

        # Materials
        "elastic": MeshEntity(ncells=10234, ncorners=4, nvertices=2361),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=146, ncorners=3, nvertices=100),
        "bc_yneg": MeshEntity(ncells=146, ncorners=3, nvertices=100),
        "bc_outer": MeshEntity(ncells=230, ncorners=3, nvertices=156),
    }


class Hex(object):
    """Mesh information for hex mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=662, ncorners=8, nvertices=1104),
        "elastic": MeshEntity(ncells=662, ncorners=8, nvertices=1104),
        "bc_xneg": MeshEntity(ncells=40, ncorners=4, nvertices=63),
        "bc_yneg": MeshEntity(ncells=40, ncorners=4, nvertices=63),
        "bc_outer": MeshEntity(ncells=64, ncorners=4, nvertices=99),
    }



# End of file
