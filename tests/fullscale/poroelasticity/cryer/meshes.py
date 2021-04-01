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
# @file tests/fullscale/poroelasticty/cryer/meshes.py
#
# @brief Mesh information for test cases.


class Tet(object):
    """
    Mesh information for tet mesh.
    """
    DOMAIN = {
        "ncells": 5398,
        "ncorners": 4,
        "nvertices": 1150,
    }
    MATERIALS = {
        "poroelastic": {
            "ncells": 5398,
            "ncorners": 4,
            "nvertices": 1150,
        }
    }
    BOUNDARIES = {
        "surface_traction": {
            "ncells": 374,
            "ncorners": 3,
            "nvertices": 212,
        },
        "surface_pressure": {
            "ncells": 374,
            "ncorners": 3,
            "nvertices": 212,
        },        
        "x_neg": {
            "ncells": 184,
            "ncorners": 3,
            "nvertices": 111,
        },
        "y_neg": {
            "ncells": 184,
            "ncorners": 3,
            "nvertices": 111,
        },
        "z_neg": {
            "ncells": 184,
            "ncorners": 3,
            "nvertices": 111,
        },
    }


class Hex(object):
    """
    Mesh information for hex mesh.
    """
    DOMAIN = {
        "ncells": 1900,
        "ncorners": 8,
        "nvertices": 2324,
    }
    MATERIALS = {
        "poroelastic": {
            "ncells": 1900,
            "ncorners": 8,
            "nvertices": 2324,
        }
    }
    BOUNDARIES = {
        "surface_traction": {
            "ncells": 300,
            "ncorners": 4,
            "nvertices": 331,
        },
        "surface_pressure": {
            "ncells": 300,
            "ncorners": 4,
            "nvertices": 331,
        },        
        "x_neg": {
            "ncells": 160,
            "ncorners": 4,
            "nvertices": 184,
        },
        "y_neg": {
            "ncells": 160,
            "ncorners": 4,
            "nvertices": 184,
        },
        "z_neg": {
            "ncells": 160,
            "ncorners": 4,
            "nvertices": 184,
        },
    }

# End of file