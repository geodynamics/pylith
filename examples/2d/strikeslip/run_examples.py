#!/usr/bin/env nemesis
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
# Copyright (c) 2010-2019 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================

from pylith.tests.FullTestApp import TestDriver, Example

import unittest

class Step01_Slip(Example):
    """
    Run 'pylith step01_slip.cfg'
    """
    NAME = "step01_slip"
    PYLITH_ARGS =  ["step01_slip.cfg"]

        
class Step02_VelBC(Example):
    """
    Run 'pylith step02_slip_velbc.cfg'.
    """
    NAME = "step02_slip_velbc"
    PYLITH_ARGS = ["step02_slip_velbc.cfg"]


class Step03_MultiSlip_VelBC(Example):
    """
    Run 'pylith step03_multislip_velbc'.
    """
    NAME = "step03_multislip_velbc"
    PYLITH_ARGS = ["step03_multislip_velbc.cfg"]

        
class ExamplesApp(TestDriver):
    """
    Driver application for running examples.
    """

    def __init__(self):
        """
        Constructor.
        """
        TestDriver.__init__(self)
        return

    def _suite(self):
        """
        Create test suite.
        """
        EXAMPLES = [
            Step01_Slip,
            Step02_VelBC,
            Step03_MultiSlip_VelBC,
            ]
        
        suite = unittest.TestSuite()

        for example in EXAMPLES:
            suite.addTest(unittest.makeSuite(example))

        return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
    Example.parse_args()
    ExamplesApp().main()


# End of file
