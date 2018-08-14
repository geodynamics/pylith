# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pyre/meshio/OutputConstraint.py
#
# @brief Python object for managing output over points with constrained degrees of freedom.
#
# Factory: output_manager

from .OutputManager import OutputManager
from .meshio import OutputConstraint as ModuleOutputConstraint


class OutputConstraint(OutputManager, ModuleOutputConstraint):
    """
    Python object for managing output over points with constrained degrees of freedom.

    Factory: observer
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputconstraint"):
        """
        Constructor.
        """
        OutputManager.__init__(self, name)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self, integrator):
        """
        Create handle to C++ object.
        """
        ModuleOutputConstraint.__init__(self, integrator)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """
    Factory associated with OutputManager.
    """
    return OutputConstraint()


# End of file
