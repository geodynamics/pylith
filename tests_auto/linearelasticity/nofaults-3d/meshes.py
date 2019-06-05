#!/usr/bin/env python
#
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
# @file tests_auto/linearelasticity/nofaults/2d/meshes.py
#
# @brief Mesh information for test cases.


class Tet(object):
    """
    Mesh information for tet mesh.
    """
    DOMAIN = {
        "ncells": 374,
        "ncorners": 4,
        "nvertices": 112,
    }
    MATERIALS = {
        "upper_crust": {
            "ncells": 140,
            "ncorners": 4,
            "nvertices": 57,
        },
        "lower_crust": {
            "ncells": 234,
            "ncorners": 4,
            "nvertices": 80,
        }
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 27,
            "ncorners": 3,
            "nvertices": 19,
        },
        "bc_xpos": {
            "ncells": 27,
            "ncorners": 3,
            "nvertices": 19,
        },
        "bc_yneg": {
            "ncells": 27,
            "ncorners": 3,
            "nvertices": 19,
        },
        "bc_ypos": {
            "ncells": 27,
            "ncorners": 3,
            "nvertices": 19,
        },
        "bc_zpos": {
            "ncells": 32,
            "ncorners": 3,
            "nvertices": 25,
        },
        "bc_zneg": {
            "ncells": 32,
            "ncorners": 3,
            "nvertices": 25,
        },
    }


class Hex(object):
    """
    Mesh information for hex mesh.
    """
    DOMAIN = {
        "ncells": 64,
        "ncorners": 8,
        "nvertices": 100,
    }
    MATERIALS = {
        "upper_crust": {
            "ncells": 16,
            "ncorners": 8,
            "nvertices": 50,
        },
        "lower_crust": {
            "ncells": 32,
            "ncorners": 8,
            "nvertices": 75,
        }
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 12,
            "ncorners": 4,
            "nvertices": 20,
        },
        "bc_xpos": {
            "ncells": 12,
            "ncorners": 4,
            "nvertices": 20,
        },
        "bc_yneg": {
            "ncells": 12,
            "ncorners": 4,
            "nvertices": 20,
        },
        "bc_ypos": {
            "ncells": 12,
            "ncorners": 4,
            "nvertices": 20,
        },
        "bc_zneg": {
            "ncells": 16,
            "ncorners": 4,
            "nvertices": 25,
        },
        "bc_zpos": {
            "ncells": 16,
            "ncorners": 4,
            "nvertices": 25,
        },
    }


# End of file
