# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.testing.FullTestApp import MeshEntity


class TetGmsh(object):
    """Mesh information for tet mesh using Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=1450, ncorners=4, nvertices=425+9),

        # Materials
        "mat_elastic": MeshEntity(ncells=1450, ncorners=4, nvertices=425+9),

        # Faults
        "fault": MeshEntity(ncells=22, ncorners=3, nvertices=18),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=84, ncorners=3, nvertices=55),
        "bc_xpos": MeshEntity(ncells=84, ncorners=3, nvertices=55),
        "bc_yneg": MeshEntity(ncells=84, ncorners=3, nvertices=55),
        "bc_ypos": MeshEntity(ncells=84, ncorners=3, nvertices=55),
        "bc_zneg": MeshEntity(ncells=162, ncorners=3, nvertices=98),
        "bc_zpos": MeshEntity(ncells=164, ncorners=3, nvertices=94),
    }


# End of file
