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

class Step01_Coseismic(Example):
    """
    Run 'pylith step01_coseismic.cfg'.
    """
    NAME = "step01_coseismic"
    PYLITH_ARGS = ["step01_coseismic.cfg"]

        
class Step02_Interseismic(Example):
    """
    Run 'pylith step02_interseismic.cfg'.
    """
    NAME = "step02_interseismic"
    PYLITH_ARGS = ["step02_interseismic.cfg"]

        
class Step03_EqCycle(Example):
    """
    Run 'pylith step03_eqcycle.cfg'.
    """
    NAME = "step03_eqcycle"
    PYLITH_ARGS = ["step03_eqcycle.cfg"]

        
class Step04_Aferslip(Example):
    """
    Run 'pylith step04_afterslip.cfg'.
    """
    NAME = "step04_afterslip"
    PYLITH_ARGS = ["step04_afterslip.cfg"]


class Step05_EqCycleSlipWeakening(Example):
    """
    Run 'pylith step05_eqcycleslipweakening.cfg'.
    """
    NAME = "step05_eqcycleslipweakening"
    PYLITH_ARGS = ["step05_eqcycleslipweakening.cfg"]


class Step06_EqCycleRateState(Example):
    """
    Run 'pylith step06_eqcycleratestate.cfg'.
    """
    NAME = "step05_eqcycleratestate"
    PYLITH_ARGS = ["step05_eqcycleratestate.cfg"]


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
            Step01_Coseismic,
            Step02_Interseismic,
            Step03_EqCycle,
            #Step04_Afterslip,
            #Step05_EqCycleSlipWeakening,
            #Step06_EqCycleRateState,
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
