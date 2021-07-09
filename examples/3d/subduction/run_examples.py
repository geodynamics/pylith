#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================

from pylith.testing.FullTestApp import TestDriver, Example

import unittest

class Step01_AxialDisp(Example):
    """Run 'pylith step01_axialdisp.cfg mat_elastic.cfg'
    """
    NAME = "step01_axialdisp"
    PYLITH_ARGS = ["step01_axialdisp.cfg", "mat_elastic.cfg"]

        
class Step02_Coseismic(Example):
    """Run 'pylith step02_coseismic.cfg'
    """
    NAME = "step02_coseismic"
    PYLITH_ARGS = ["step02_coseismic.cfg", "mat_viscoelastic.cfg", "solver_fieldsplit.cfg"]

        
class Step03_Interseismic(Example):
    """Run 'pylith step03_interseismic.cfg mat_viscoelastic.cfg solver_fieldsplit.cfg'
    """
    NAME = "step03_interseismic"
    PYLITH_ARGS = ["step03_interseismic.cfg", "mat_viscoelastic.cfg", "solver_fieldsplit.cfg"]

        
class Step04_EqCycle(Example):
    """Run 'pylith step04_eqcycle.cfg'.
    """
    NAME = "step04_eqcycle"
    PYLITH_ARGS = ["step04_eqcycle.cfg", "mat_viscoelastic.cfg", "solver_fieldsplit.cfg"]

        
class Step05_EqCycle_SlipWeakening(Example):
    """Run 'pylith step05_eqcycle_slipweakening.cfg mat_elastic.cfg solver_fieldsplit.cfg'
    """
    NAME = "step05_eqcycle_slipweakening"
    PYLITH_ARGS = ["step05_eqcycle_slipweakening.cfg", "mat_elastic.cfg", "solver_fieldsplit.cfg"]

        
class Step06_SlowSlip(Example):
    """Run 'pylith step06_slowslip.cfg mat_elastic.cfg solver_fieldsplit.cfg'
    """
    NAME = "step06_slowslip"
    PYLITH_ARGS = ["step06_slowslip.cfg", "mat_elastic.cfg", "solver_fieldsplit.cfg"]

        
class Step07a_GreensFns_LeftLateral(Example):
    """Run 'pylith step07a_leftlateral.cfg mat_elastic.cfg solver_fieldsplit.cfg'
    """
    NAME = "step07a_greensfns_leftlateral"
    PYLITH_ARGS = [
        "--problem=pylith.problems.GreensFns",
        "step07a_leftlateral.cfg",
        "mat_elastic.cfg",
        "solver_fieldsplit.cfg"
        ]

        
class Step07b_GreensFns_Reverse(Example):
    """Run 'pylith step07b_reverse.cfg mat_elastic.cfg solver_fieldsplit.cfg'
    """
    NAME = "step07b_greensfns_reverse"
    PYLITH_ARGS = [
        "--problem=pylith.problems.GreensFns",
        "step07b_reverse.cfg",
        "mat_elastic.cfg",
        "solver_fieldsplit.cfg"
        ]

        
class Step08a_Gravity_RefState(Example):
    """Run 'pylith step08a_gravity_refstate.cfg mat_elastic.cfg solver_algebraicmultigrid.cfg'
    """
    NAME = "step08a_gravity_refstate"
    PYLITH_ARGS = ["step08a_gravity_refstate.cfg", "mat_elastic.cfg", "solver_algebraicmultigrid.cfg"]

        
class Step08b_Gravity_Incompressible(Example):
    """Run 'pylith step08b_gravity_incompressible.cfg mat_incompressible.cfg solver_algebraicmultigrid.cfg'
    """
    NAME = "step08b_gravity_incompressible"
    PYLITH_ARGS = ["step08b_gravity_incompressible.cfg", "mat_incompressible.cfg", "solver_algebraicmultigrid.cfg"]

        
class Step08c_Gravity_Viscoelastic(Example):
    """Run 'pylith step08c_gravity_viscoelastic.cfg mat_viscoelastic.cfg solver_algebraicmultigrid.cfg'
    """
    NAME = "step08c_gravity_viscoelastic"
    PYLITH_ARGS = ["step08c_gravity_viscoelastic.cfg", "mat_viscoelastic.cfg", "solver_algebraicmultigrid.cfg"]


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
            Step02_Coseismic,
            Step03_Interseismic,
            Step04_EqCycle,
            #Step05_EqCycle_SlipWeakening,
            Step06_SlowSlip,
            #Step07a_GreensFns_LeftLateral,
            #Step07b_GreensFns_Reverse,
            Step08a_Gravity_RefState,
            Step08b_Gravity_Incompressible,
            Step08c_Gravity_Viscoelastic,
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
