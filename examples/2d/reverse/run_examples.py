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

class Step01_Gravity(Example):
    """
    Run 'pylith step01_gravity.cfg'.
    """
    NAME = "step01_gravity"
    PYLITH_ARGS = ["step01_gravity.cfg"]

        
class Step02_Gravity_RefState(Example):
    """
    Run 'pylith step02_gravity_refstate.cfg'.
    """
    NAME = "step02_gravity_refstate"
    PYLITH_ARGS = ["step02_gravity_refstate.cfg"]

        
class Step03_Gravity_Incompressible(Example):
    """
    Run 'pylith step03_gravity_incompressible.cfg'.
    """
    NAME = "step03_gravity_incompressible"
    PYLITH_ARGS = ["step03_gravity_incompressible.cfg"]

        
class Step04_SurfLoad(Example):
    """
    Run 'pylith step04_surfload.cfg'.
    """
    NAME = "step04_surfload"
    PYLITH_ARGS = ["step04_surfload.cfg"]


class Step05_OneFault(Example):
    """
    Run 'pylith step05_onefault.cfg solver_onefault.cfg'.
    """
    NAME = "step05_onefault"
    PYLITH_ARGS = ["step05_onefault.cfg", "solver_onefault.cfg"]


class Step06_TwoFaults_Elastic(Example):
    """
    Run 'pylith step06_twofaults.cfg solver_twofaults.cfg'.
    """
    NAME = "step06_twofaults_elastic"
    PYLITH_ARGS = ["step06_twofaults_elastic.cfg", "solver_twofaults.cfg"]


class Step07_TwoFaults_Maxwell(Example):
    """
    Run 'pylith step07_twofaults_maxwell.cfg solver_twofaults.cfg'.
    """
    NAME = "step07_twofaults_maxwell"
    PYLITH_ARGS = ["step07_twofaults_maxwell.cfg", "solver_twofaults.cfg"]


class Step08_TwoFaults_Powerlaw(Example):
    """
    Run 'pylith step08_twofaults_powerlaw.cfg solver_twofaults.cfg'.
    """
    NAME = "step08_twofaults_powerlaw"
    PYLITH_ARGS = ["step08_twofaults_powerlaw.cfg", "solver_twofaults.cfg"]


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
            Step01_Gravity,
            Step02_Gravity_RefState,
            Step03_Gravity_Incompressible,
            Step04_SurfLoad,
            Step05_OneFault,
            Step06_TwoFaults_Elastic,
            #Step07_TwoFaults_Maxwell,
            #Step08_TwoFaults_Powerlaw,
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
