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
# @file pythia.pyre/meshio/OutputSoln.py
#
# @brief Python object for managing output of solution information over the domain.
#
# Factory: observer

from .OutputSoln import OutputSoln
from .meshio import OutputSolnDomain as ModuleOutputSolnDomain


class OutputSolnDomain(OutputSoln, ModuleOutputSolnDomain):
    """Python object for managing output of finite-element solution
    information.

    FACTORY: observer
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputsolndomain"):
        """Constructor.
        """
        OutputSoln.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Do mimimal initialization.
        """
        OutputSoln.preinitialize(self, problem)

        identifier = self.aliases[-1]
        self.writer.setFilename(problem.defaults.outputDir, problem.defaults.simName, identifier)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        OutputSoln._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleOutputSolnDomain.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """Factory associated with OutputSolnDomain.
    """
    return OutputSolnDomain()


# End of file
