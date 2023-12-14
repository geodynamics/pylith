# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pythia.pyre.components.Component import Component


def validateName(value):
    if not value.strip():
        raise ValueError("Missing required property 'name' in default options for problem.")
    return value.strip()


class ProblemDefaults(Component):
    """
    Default options for a problem.
    Specifying defaults at the problem level (here) will override defaults for individual components.
    Non-default values specified for individual components will override the problem defaults (specified here).
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.defaults]
            output_directory = output
            name = step01
            quadrature_order = 1
            output_basis_order = 0
        """
    }

    import pythia.pyre.inventory

    outputDir = pythia.pyre.inventory.str("output_directory", default="output")
    outputDir.meta['tip'] = "Directory for output."

    simName = pythia.pyre.inventory.str("name", default="", validator=validateName)
    simName.meta['tip'] = "Name for the problem (used with output_directory for default output filenames)."

    quadOrder = pythia.pyre.inventory.int("quadrature_order", default=1, validator=pythia.pyre.inventory.greater(0))
    quadOrder.meta['tip'] = "Finite-element quadrature order."

    outputBasisOrder = pythia.pyre.inventory.int("output_basis_order", default=1, validator=pythia.pyre.inventory.choice([0,1]))
    outputBasisOrder.meta['tip'] = "Default basis order for output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="problem_defaults"):
        """Constructor.
        """
        Component.__init__(self, name, facility="problem_defaults")
        return

    def preinitialize(self):
        """Do minimal initialization.
        """
        return


# FACTORIES ////////////////////////////////////////////////////////////


def problem_defaults():
    """Factory associated with ProblemDefaults.
    """
    return ProblemDefaults()


# End of file
