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


class CylinderTet(object):
    """
    Mesh information for cylindrical tet mesh.
    """
    DOMAIN = {
        "ncells": 84,
        "ncorners": 4,
        "nvertices": 38,
    }
    MATERIALS = {
        "elastic": {
            "ncells": 84,
            "ncorners": 4,
            "nvertices": 38,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 7,
            "ncorners": 3,
            "nvertices": 8,
        },
        "bc_yneg": {
            "ncells": 7,
            "ncorners": 3,
            "nvertices": 8,
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
            "ncells": 8,
            "ncorners": 3,
            "nvertices": 10,
        },
        "bc_inner": {
            "ncells": 10,
            "ncorners": 3,
            "nvertices": 11,
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


class CylinderHex(object):
    """
    Mesh information for cylindrical hex mesh.
    """
    DOMAIN = {
        "ncells": 6,
        "ncorners": 8,
        "nvertices": 24,
    }
    MATERIALS = {
        "elastic": {
            "ncells": 6,
            "ncorners": 8,
            "nvertices": 24,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 2,
            "ncorners": 4,
            "nvertices": 6,
        },
        "bc_yneg": {
            "ncells": 2,
            "ncorners": 4,
            "nvertices": 6,
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
            "ncells": 3,
            "ncorners": 4,
            "nvertices": 8,
        },
        "bc_inner": {
            "ncells": 3,
            "ncorners": 4,
            "nvertices": 8,
        },
    }


# End of file
