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
# @file tests/fullscale/linearelasticity/cylinder-2d/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity

class Tri(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=757, ncorners=3, nvertices=415),

        # Materials
        "elastic": MeshEntity(ncells=757, ncorners=3, nvertices=415),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "bc_yneg": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "bc_outer": MeshEntity(ncells=31, ncorners=2, nvertices=32),
    }


class Quad(object):
    """Mesh information for quad mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=331, ncorners=4, nvertices=368),
        "elastic": MeshEntity(ncells=331, ncorners=4, nvertices=368),
        "bc_xneg": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "bc_yneg": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "bc_outer": MeshEntity(ncells=32, ncorners=2, nvertices=33),
    }



# End of file
