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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

import numpy

def test_double(obj, valuesE, values):
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
            obj.assertAlmostEqual(v, vE, 6)
        else:
            obj.assertAlmostEqual(1.0, v/vE, 6)

    return


# End of file
