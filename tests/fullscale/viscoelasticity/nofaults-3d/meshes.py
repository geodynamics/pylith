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
# @file tests/fullscale/linearelasticity/nofaults/3d/meshes.py
#
# @brief Mesh information for test cases.


class Tet(object):
    """
    Mesh information for tri mesh.
    """
    DOMAIN = {
        "ncells": 530,
        "ncorners": 4,
        "nvertices": 149,
    }
    MATERIALS = {
        "viscomat": {
            "ncells": 530,
            "ncorners": 4,
            "nvertices": 149,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 32,
            "ncorners": 3,
            "nvertices": 25,
        },
        "bc_xpos": {
            "ncells": 32,
            "ncorners": 3,
            "nvertices": 25,
        },
        "bc_yneg": {
            "ncells": 32,
            "ncorners": 3,
            "nvertices": 25,
        },
        "bc_ypos": {
            "ncells": 32,
            "ncorners": 3,
            "nvertices": 25,
        },
        "bc_zneg": {
            "ncells": 32,
            "ncorners": 3,
            "nvertices": 25,
        },
        "bc_zpos": {
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
        "nvertices": 125,
    }
    MATERIALS = {
        "viscomat": {
            "ncells": 64,
            "ncorners": 8,
            "nvertices": 125,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 16,
            "ncorners": 4,
            "nvertices": 25,
        },
        "bc_xpos": {
            "ncells": 16,
            "ncorners": 4,
            "nvertices": 25,
        },
        "bc_yneg": {
            "ncells": 16,
            "ncorners": 4,
            "nvertices": 25,
        },
        "bc_ypos": {
            "ncells": 16,
            "ncorners": 4,
            "nvertices": 25,
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
