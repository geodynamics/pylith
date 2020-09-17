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
# @file tests/fullscale/cornercases/3d/meshes.py
#
# @brief Mesh information for test cases.


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


# End of file
