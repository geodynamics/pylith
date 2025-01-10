# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
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
