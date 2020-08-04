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
# @file tests/fullscale/poroelasticty/terzaghi/meshes.py
#
# @brief Mesh information for test cases.


class Tri(object):
    """
    Mesh information for tri mesh.
    """
    DOMAIN = {
        "ncells": 266,
        "ncorners": 3,
        "nvertices": 170,
    }
    MATERIALS = {
        "poroelastic": {
            "ncells": 266,
            "ncorners": 3,
            "nvertices": 170,
        }
    }
    BOUNDARIES = {
        "x_neg": {
            "ncells": 3,
            "ncorners": 2,
            "nvertices": 4,
        },
        "x_pos": {
            "ncells": 3,
            "ncorners": 2,
            "nvertices": 4,
        },
        "y_neg": {
            "ncells": 33,
            "ncorners": 2,
            "nvertices": 34,
        },
        "y_pos": {
            "ncells": 33,
            "ncorners": 2,
            "nvertices": 34,
        },
    }


class Quad(object):
    """
    Mesh information for quad mesh.
    """
    DOMAIN = {
        "ncells": 160,
        "ncorners": 4,
        "nvertices": 205,
    }
    MATERIALS = {
        "poroelastic": {
            "ncells": 100,
            "ncorners": 4,
            "nvertices": 121,
        }
    }
    BOUNDARIES = {
        "x_neg": {
            "ncells": 4,
            "ncorners": 2,
            "nvertices": 5,
        },
        "x_pos": {
            "ncells": 4,
            "ncorners": 2,
            "nvertices": 5,
        },
        "y_neg": {
            "ncells": 40,
            "ncorners": 2,
            "nvertices": 41,
        },
        "y_pos": {
            "ncells": 40,
            "ncorners": 2,
            "nvertices": 41,
        },
    }


# End of file
