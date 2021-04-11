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

from pylith.testing.FullTestApp import TestDriver, Example

import unittest

class Step01_AxialDisp(Example):
    """Run 'pylith step01_axialdisp.cfg'
    """
    NAME = "step01_axialdisp"
    PYLITH_ARGS = ["step01_axialdisp.cfg"]

        
class Step02_ShearDisp(Example):
    """Run 'pylith step02_sheardisp.cfg'
    """
    NAME = "step02_sheardisp"
    PYLITH_ARGS = ["step02_sheardisp.cfg"]

        
class Step03_ShearDispTract(Example):
    """Run 'pylith step03_sheardisptract.cfg'
    """
    NAME = "step03_sheardisptract"
    PYLITH_ARGS = ["step03_sheardisptract.cfg"]

        
class Step04_ShearDispIC(Example):
    """Run 'pylith step04_sheardispic.cfg'.
    """
    NAME = "step04_sheardispic"
    PYLITH_ARGS = ["step04_sheardispic.cfg"]

        
class Step05_ShearDispTractRate(Example):
    """Run 'pylith step05_sheardisptractrate.cfg'
    """
    NAME = "step05_sheardisptractrate"
    PYLITH_ARGS = ["step05_sheardisptractrate.cfg"]

        
class ExamplesApp(TestDriver):
    """Driver application for running examples.
    """

    def __init__(self):
        """Constructor.
        """
        TestDriver.__init__(self)
        return

    def _suite(self):
        """Create test suite.
        """
        EXAMPLES = [
            Step01_AxialDisp,
            Step02_ShearDisp,
            Step03_ShearDispTract,
            Step04_ShearDispIC,
            Step05_ShearDispTractRate,
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
