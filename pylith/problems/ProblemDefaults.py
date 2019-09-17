# ----------------------------------------------------------------------
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
# ----------------------------------------------------------------------
#
# @file pylith/problems/ProblemDefaults.py
#
# @brief Python class for default options for a problem.
#
# Factory: problem_defaults.

from pyre.components.Component import Component


def validateName(value):
    if not value.strip():
        raise ValueError("Missing problem name.")
    return value.strip()


class ProblemDefaults(Component):
    """
    Python class for default options for a problem.

    INVENTORY

    Properties
      - *output_directory* Directory for output.
      - *name* Name of problem.
      - *quadrature_order* Finite-element quadrature order.

    Facilities
      None
    """

    import pyre.inventory

    outputDir = pyre.inventory.str("output_directory", default="output")
    outputDir.meta['tip'] = "Directory for output."

    name = pyre.inventory.str("name", default="", validator=validateName)
    name.meta['tip'] = "Name for the problem (used with output_directory for default output filenames)."

    quadOrder = pyre.inventory.int("quadrature_order", default=1, validator=pyre.inventory.greater(0))
    quadOrder.meta['tip'] = "Finite-element quadrature order."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="problem_defaults"):
        """
        Constructor.
        """
        Component.__init__(self, name, facility="problem_defaults")
        return

    def preinitialize(self, mesh):
        """
        Do minimal initialization.
        """
        return


# FACTORIES ////////////////////////////////////////////////////////////


def problem_defaults():
    """
    Factory associated with ProblemDefaults.
    """
    return ProblemDefaults()


# End of file
