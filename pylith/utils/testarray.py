#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
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
        obj.assertAlmostEqual(1.0, v/vE, 6)

    return


# End of file
