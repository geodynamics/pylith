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
        "domain": MeshEntity(ncells=124, ncorners=3, nvertices=79),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "elastic_xpos": MeshEntity(ncells=64, ncorners=3, nvertices=45),
        "elastic_xneg": MeshEntity(ncells=60, ncorners=3, nvertices=43),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
    }


class Quad(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=64, ncorners=4, nvertices=81),
        "elastic_xpos": MeshEntity(ncells=32, ncorners=4, nvertices=45),
        "elastic_xneg": MeshEntity(ncells=32, ncorners=4, nvertices=45),
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),
    }



# End of file
