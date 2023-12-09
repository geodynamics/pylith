# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/fullscale/poroelasticty/terzaghi/meshes.py
#
# @brief Mesh information for test cases.

from pylith.testing.FullTestApp import MeshEntity


class Tri(object):
    """Mesh information for tri mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=902, ncorners=3, nvertices=492),

        # Materials
        "poroelastic": MeshEntity(ncells=902, ncorners=3, nvertices=492),

        # Boundaries
        "x_neg": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "x_pos": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "y_neg": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "y_pos_dir": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "y_pos_neu": MeshEntity(ncells=20, ncorners=2, nvertices=21),
    }


class Quad(object):
    """Mesh information for quad mesh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=400, ncorners=4, nvertices=441),

        # Materials
        "poroelastic": MeshEntity(ncells=400, ncorners=4, nvertices=441),

        # Boundaries
        "x_neg": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "x_pos": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "y_neg": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "y_pos_dir": MeshEntity(ncells=20, ncorners=2, nvertices=21),
        "y_pos_neu": MeshEntity(ncells=20, ncorners=2, nvertices=21),
    }


# End of file
