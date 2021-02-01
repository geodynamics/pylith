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
        "ncells": 902,
        "ncorners": 3,
        "nvertices": 492,
    }
    MATERIALS = {
        "poroelastic": {
            "ncells": 902,
            "ncorners": 3,
            "nvertices": 492,
        }
    }
    BOUNDARIES = {
        "x_neg": {
            "ncells": 20,
            "ncorners": 2,
            "nvertices": 21,
        },
        "x_pos": {
            "ncells": 20,
            "ncorners": 2,
            "nvertices": 21,
        },
        "y_pos_neu": {
            "ncells": 20,
            "ncorners": 2,
            "nvertices": 21,
        },
        "y_pos_dir": {
            "ncells": 20,
            "ncorners": 2,
            "nvertices": 21,
        },
        "y_neg": {
            "ncells": 20,
            "ncorners": 2,
            "nvertices": 21,
        },
    }


class Quad(object):
    """
    Mesh information for quad mesh.
    """
    DOMAIN = {
        "ncells": 400,
        "ncorners": 4,
        "nvertices": 441,
    }
    MATERIALS = {
        "poroelastic": {
            "ncells": 400,
            "ncorners": 4,
            "nvertices": 441,
        }
    }
    BOUNDARIES = {
        "x_neg": {
            "ncells": 20,
            "ncorners": 2,
            "nvertices": 21,
        },
        "x_pos": {
            "ncells": 20,
            "ncorners": 2,
            "nvertices": 21,
        },
        "y_pos_neu": {
            "ncells": 20,
            "ncorners": 2,
            "nvertices": 21,
        },
        "y_pos_dir": {
            "ncells": 20,
            "ncorners": 2,
            "nvertices": 21,
        },
        "y_neg": {
            "ncells": 20,
            "ncorners": 2,
            "nvertices": 21,
        },
    }


# End of file
