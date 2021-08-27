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
#
# @file tests/fullscale/linearelasticity/nofaults/2d/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity


class Tri(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=124, ncorners=3, nvertices=88),

        # Materials
        "mat_xneg": MeshEntity(ncells=30, ncorners=3, nvertices=26),
        "mat_xmid": MeshEntity(ncells=30, ncorners=3, nvertices=26),
        "mat_xposypos": MeshEntity(ncells=32, ncorners=3, nvertices=25),
        "mat_xposyneg": MeshEntity(ncells=32, ncorners=3, nvertices=25),

        # Faults
        "fault": MeshEntity(ncells=8, ncorners=2, nvertices=9),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=10),
    }


class Quad(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=64, ncorners=4, nvertices=90),

        # Materials
        "mat_xneg": MeshEntity(ncells=16, ncorners=3, nvertices=27),
        "mat_xmid": MeshEntity(ncells=16, ncorners=3, nvertices=27),
        "mat_xposypos": MeshEntity(ncells=16, ncorners=3, nvertices=25),
        "mat_xposyneg": MeshEntity(ncells=16, ncorners=3, nvertices=25),

        # Faults
        "fault": MeshEntity(ncells=8, ncorners=2, nvertices=9),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=10),
    }


# End of file
