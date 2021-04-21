# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/cornercases/3d/meshes.py
#
# @brief Mesh information for test cases.


class TetPgram(object):
    """
    Mesh information for cylindrical tet mesh.
    """
    DOMAIN = {
        "ncells": 96,
        "ncorners": 4,
        "nvertices": 35,
    }
    MATERIALS = {
        "elastic": {
            "ncells": 96,
            "ncorners": 4,
            "nvertices": 35,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 8,
            "ncorners": 3,
            "nvertices": 9,
        },
        "bc_xpos": {
            "ncells": 8,
            "ncorners": 3,
            "nvertices": 9,
        },
        "bc_yneg": {
            "ncells": 8,
            "ncorners": 3,
            "nvertices": 9,
        },
        "bc_ypos": {
            "ncells": 8,
            "ncorners": 3,
            "nvertices": 9,
        },
        "bc_zneg": {
            "ncells": 8,
            "ncorners": 3,
            "nvertices": 9,
        },
        "bc_zpos": {
            "ncells": 8,
            "ncorners": 3,
            "nvertices": 9,
        },
        "bc_outer": {
            "ncells": 48,
            "ncorners": 3,
            "nvertices": 26,
        },
    }


class TetCylinder(object):
    """
    Mesh information for cylindrical tet mesh.
    """
    DOMAIN = {
        "ncells": 139,
        "ncorners": 4,
        "nvertices": 53,
    }
    MATERIALS = {
        "elastic": {
            "ncells": 139,
            "ncorners": 4,
            "nvertices": 53,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 9,
            "ncorners": 3,
            "nvertices": 10,
        },
        "bc_yneg": {
            "ncells": 9,
            "ncorners": 3,
            "nvertices": 10,
        },
        "bc_zneg": {
            "ncells": 17,
            "ncorners": 3,
            "nvertices": 15,
        },
        "bc_zpos": {
            "ncells": 15,
            "ncorners": 3,
            "nvertices": 14,
        },
        "bc_outer": {
            "ncells": 14,
            "ncorners": 3,
            "nvertices": 14,
        },
        "bc_inner": {
            "ncells": 18,
            "ncorners": 3,
            "nvertices": 16,
        },
    }


class Tet(object):
    """
    Mesh information for tet mesh.
    """
    DOMAIN = {
        "ncells": 5,
        "ncorners": 4,
        "nvertices": 8,
    }
    MATERIALS = {
        "elastic": {
            "ncells": 5,
            "ncorners": 4,
            "nvertices": 8,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 2,
            "ncorners": 3,
            "nvertices": 4,
        },
        "bc_xpos": {
            "ncells": 2,
            "ncorners": 3,
            "nvertices": 4,
        },
    }


class HexPgram(object):
    """
    Mesh information for parallelogram hex mesh.
    """
    DOMAIN = {
        "ncells": 8,
        "ncorners": 8,
        "nvertices": 27,
    }
    MATERIALS = {
        "elastic": {
            "ncells": 8,
            "ncorners": 8,
            "nvertices": 27,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 4,
            "ncorners": 4,
            "nvertices": 9,
        },
        "bc_xpos": {
            "ncells": 4,
            "ncorners": 4,
            "nvertices": 9,
        },
        "bc_yneg": {
            "ncells": 4,
            "ncorners": 4,
            "nvertices": 9,
        },
        "bc_ypos": {
            "ncells": 4,
            "ncorners": 4,
            "nvertices": 9,
        },
        "bc_zneg": {
            "ncells": 4,
            "ncorners": 4,
            "nvertices": 9,
        },
        "bc_zpos": {
            "ncells": 4,
            "ncorners": 4,
            "nvertices": 9,
        },
        "bc_outer": {
            "ncells": 24,
            "ncorners": 4,
            "nvertices": 26,
        },
    }


class Hex(object):
    """
    Mesh information for hex mesh.
    """
    DOMAIN = {
        "ncells": 1,
        "ncorners": 8,
        "nvertices": 8,
    }
    MATERIALS = {
        "elastic": {
            "ncells": 1,
            "ncorners": 8,
            "nvertices": 8,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 1,
            "ncorners": 4,
            "nvertices": 4,
        },
        "bc_xpos": {
            "ncells": 1,
            "ncorners": 4,
            "nvertices": 4,
        },
        "bc_yneg": {
            "ncells": 1,
            "ncorners": 4,
            "nvertices": 4,
        },
        "bc_ypos": {
            "ncells": 1,
            "ncorners": 4,
            "nvertices": 4,
        },
        "bc_zneg": {
            "ncells": 1,
            "ncorners": 4,
            "nvertices": 4,
        },
        "bc_zpos": {
            "ncells": 1,
            "ncorners": 4,
            "nvertices": 4,
        },
        "bc_domain": {
            "ncells": 6,
            "ncorners": 4,
            "nvertices": 8,
        },
    }


class HexCylinder(object):
    """
    Mesh information for cylindrical hex mesh.
    """
    DOMAIN = {
        "ncells": 12,
        "ncorners": 8,
        "nvertices": 36,
    }
    MATERIALS = {
        "elastic": {
            "ncells": 12,
            "ncorners": 8,
            "nvertices": 36,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 4,
            "ncorners": 4,
            "nvertices": 9,
        },
        "bc_yneg": {
            "ncells": 4,
            "ncorners": 4,
            "nvertices": 9,
        },
        "bc_zneg": {
            "ncells": 6,
            "ncorners": 4,
            "nvertices": 12,
        },
        "bc_zpos": {
            "ncells": 6,
            "ncorners": 4,
            "nvertices": 12,
        },
        "bc_outer": {
            "ncells": 6,
            "ncorners": 4,
            "nvertices": 12,
        },
        "bc_inner": {
            "ncells": 6,
            "ncorners": 4,
            "nvertices": 12,
        },
    }


# End of file
