# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.testing.FullTestApp import MeshEntity


class TriGmsh(object):
    """Mesh information for tri mesh using Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=164, ncorners=3, nvertices=99+2*9),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "mat_xneg": MeshEntity(ncells=38, ncorners=3, nvertices=30),
        "mat_xmid": MeshEntity(ncells=38, ncorners=3, nvertices=30),
        "mat_xposypos": MeshEntity(ncells=44, ncorners=3, nvertices=31),
        "mat_xposyneg": MeshEntity(ncells=44, ncorners=3, nvertices=31),

        # Faults
        "fault_xmid": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "fault_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9+2),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9+2),
        "boundary_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9+2),
    }


class QuadGmsh(object):
    """Mesh information for quad mesh using Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=100, ncorners=4, nvertices=121+2*11),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "mat_xneg": MeshEntity(ncells=30, ncorners=3, nvertices=44),
        "mat_xmid": MeshEntity(ncells=20, ncorners=3, nvertices=33),
        "mat_xposypos": MeshEntity(ncells=25, ncorners=3, nvertices=36),
        "mat_xposyneg": MeshEntity(ncells=25, ncorners=3, nvertices=36),

        # Faults
        "fault_xmid": MeshEntity(ncells=10, ncorners=2, nvertices=11),
        "fault_xneg": MeshEntity(ncells=10, ncorners=2, nvertices=11),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=10, ncorners=2, nvertices=11),
        "bc_xpos": MeshEntity(ncells=10, ncorners=2, nvertices=11),
        "bc_yneg": MeshEntity(ncells=11, ncorners=2, nvertices=11+2),
        "bc_ypos": MeshEntity(ncells=11, ncorners=2, nvertices=11+2),
        "boundary_ypos": MeshEntity(ncells=11, ncorners=2, nvertices=11+2),
    }


class TriCubit(object):
    """Mesh information for tri mesh using Cubit.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=124, ncorners=3, nvertices=79+2*9),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "mat_xneg": MeshEntity(ncells=30, ncorners=3, nvertices=26),
        "mat_xmid": MeshEntity(ncells=30, ncorners=3, nvertices=26),
        "mat_xposypos": MeshEntity(ncells=32, ncorners=3, nvertices=25),
        "mat_xposyneg": MeshEntity(ncells=32, ncorners=3, nvertices=25),

        # Faults
        "fault_xmid": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "fault_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9+2),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9+2),
        "boundary_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9+2),
    }


class QuadCubit(object):
    """Mesh information for quad mesh using Cubit.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=64, ncorners=4, nvertices=81+2*9),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "mat_xneg": MeshEntity(ncells=16, ncorners=3, nvertices=27),
        "mat_xmid": MeshEntity(ncells=16, ncorners=3, nvertices=27),
        "mat_xposypos": MeshEntity(ncells=16, ncorners=3, nvertices=25),
        "mat_xposyneg": MeshEntity(ncells=16, ncorners=3, nvertices=25),

        # Faults
        "fault_xmid": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "fault_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9+2),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9+2),
        "boundary_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9+2),
    }


# End of file
