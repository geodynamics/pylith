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
# @file tests/fullscale/linearelasticity/nofaults/3d/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity


class Tet(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=530, ncorners=4, nvertices=149),

        # Materials
        "viscomat": MeshEntity(ncells=530, ncorners=4, nvertices=149),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=32, ncorners=3, nvertices=25),
        "bc_xpos": MeshEntity(ncells=32, ncorners=3, nvertices=25),
        "bc_yneg": MeshEntity(ncells=32, ncorners=3, nvertices=25),
        "bc_ypos": MeshEntity(ncells=32, ncorners=3, nvertices=25),
        "bc_zneg": MeshEntity(ncells=32, ncorners=3, nvertices=25),
        "bc_zpos": MeshEntity(ncells=32, ncorners=3, nvertices=25),
    }


class Hex(object):
    """Mesh information for hex mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=64, ncorners=8, nvertices=125),

        # Materials
        "viscomat": MeshEntity(ncells=64, ncorners=8, nvertices=125),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=16, ncorners=4, nvertices=25),
        "bc_xpos": MeshEntity(ncells=16, ncorners=4, nvertices=25),
        "bc_yneg": MeshEntity(ncells=16, ncorners=4, nvertices=25),
        "bc_ypos": MeshEntity(ncells=16, ncorners=4, nvertices=25),
        "bc_zneg": MeshEntity(ncells=16, ncorners=4, nvertices=25),
        "bc_zpos": MeshEntity(ncells=16, ncorners=4, nvertices=25),
    }


# End of file
