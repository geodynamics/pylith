#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

import numpy

def test_scalararray(obj, valuesE, values, places=6):
    """
    Check whether arrays containing double values match.

    @param obj Test object
    @valuesE Array of expected values
    @values Array of values to check
    """

    # Check shape
    obj.assertEqual(len(valuesE.shape), len(values.shape))
    for (dimE, dim) in zip(valuesE.shape, values.shape):
        obj.assertEqual(dimE, dim)

    # Check values
    for (vE, v) in zip(numpy.reshape(valuesE, -1),
                       numpy.reshape(values, -1)):
        if vE == 0.0:
            obj.assertAlmostEqual(v, vE, places)
        else:
            obj.assertAlmostEqual(1.0, v/vE, places)

    return


# End of file
