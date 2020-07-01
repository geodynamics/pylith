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
# @file tests/fullscale/linearelasticity/nofaults/2d/meshes.py
#
# @brief Mesh information for test cases.


class Tri(object):
    """
    Mesh information for tri mesh.
    """
    DOMAIN = {
        "ncells": 150,
        "ncorners": 3,
        "nvertices": 92,
    }
    MATERIALS = {
        "viscomat": {
            "ncells": 150,
            "ncorners": 3,
            "nvertices": 92,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_xpos": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_yneg": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_ypos": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
    }


class Quad(object):
    """
    Mesh information for tri mesh.
    """
    DOMAIN = {
        "ncells": 64,
        "ncorners": 4,
        "nvertices": 81,
    }
    MATERIALS = {
        "viscomat": {
            "ncells": 64,
            "ncorners": 4,
            "nvertices": 81,
        },
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_xpos": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_yneg": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_ypos": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
    }


# End of file
