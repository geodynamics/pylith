# ----------------------------------------------------------------------
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
# ----------------------------------------------------------------------
#
# @file pylith/problems/ProblemDefaults.py
#
# @brief Python class for default options for a problem.
#
# Factory: problem_defaults.

from pythia.pyre.components.Component import Component


def validateName(value):
    if not value.strip():
        raise ValueError("Missing required property 'name' in default options for problem.")
    return value.strip()


class ProblemDefaults(Component):
    """Python class for default options for a problem.

    FACTORY: problem_defaults
    """

    import pythia.pyre.inventory

    outputDir = pythia.pyre.inventory.str("output_directory", default="output")
    outputDir.meta['tip'] = "Directory for output."

    simName = pythia.pyre.inventory.str("name", default="", validator=validateName)
    simName.meta['tip'] = "Name for the problem (used with output_directory for default output filenames)."

    quadOrder = pythia.pyre.inventory.int("quadrature_order", default=1, validator=pythia.pyre.inventory.greater(0))
    quadOrder.meta['tip'] = "Finite-element quadrature order."

    outputBasisOrder = pythia.pyre.inventory.int("output_basis_order", default=1, validator=pythia.pyre.inventory.choice((0,1)))
    outputBasisOrder.meta['tip'] = "Basis order for output."

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
