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
##
# @file tests_auto/linearelasticity/nofaults/TestTri.py
##
# @brief Base class for tests with 2-D simplex mesh.

import unittest


class TestQuad(unittest.TestCase):
    """
    Generic tests for problems using 2-D mesh.
    """

    def setUp(self):
        """
        Setup for tests.
        """
        self.domain = {
            "ncells": 64,
            "ncorners": 4,
            "nvertices": 81,
        }
        self.materials = {
            "elastic_xpos": {
                "ncells": 32,
                "ncorners": 4,
                "nvertices": 36,
            },
            "elastic_xneg": {
                "ncells": 32,
                "ncorners": 4,
                "nvertices": 36,
            }
        }


# End of file
