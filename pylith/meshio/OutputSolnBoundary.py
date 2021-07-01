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
# @file pythia.pyre/meshio/OutputSolnBoundary.py
#
# @brief Python object for managing output of finite-element solution
# information over a subdomain.
#
# Factory: observer

from .OutputSoln import OutputSoln
from .meshio import OutputSolnBoundary as ModuleOutputSolnBoundary


def validateLabel(value):
    """Validate label for group/nodeset/pset.
    """
    if not value.strip():
        raise ValueError("Label for group/nodeset/pset in mesh not specified.")
    return value


class OutputSolnBoundary(OutputSoln, ModuleOutputSolnBoundary):
    """Python object for managing output of finite-element solution
    information over a boundary.

    Factory: observer
    """

    import pythia.pyre.inventory

    label = pythia.pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for boundary."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputsolnsubset"):
        """Constructor.
        """
        OutputSoln.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Do mimimal initialization.
        """
        OutputSoln.preinitialize(self, problem)
        ModuleOutputSolnBoundary.setLabel(self, self.label)

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
        ModuleOutputSolnBoundary.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """Factory associated with OutputSoln.
    """
    return OutputSolnBoundary()


# End of file
