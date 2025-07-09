# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .OutputSoln import OutputSoln
from .meshio import OutputSolnDomain as ModuleOutputSolnDomain


class OutputSolnDomain(OutputSoln, ModuleOutputSolnDomain):
    """
    Output of solution subfields over the simulation domain.

    :::{tip}
    Most output information can be configured at the problem level using the [`ProblemDefaults` Component](../problems/ProblemDefaults.md).
    :::

    Implements `OutputSoln`.
    """
    DOC_CONFIG = {
        "cfg": """
            [observer]
            data_fields = [displacement]

            # Skip two time steps between output.
            trigger = pylith.meshio.OutputTriggerStep
            trigger.num_skip = 2

            # Write output to HDF5 file with name `domain.h5`.
            writer = pylith.meshio.DataWriterHDF5
            writer.filename = domain.h5

            # Output with a basis order of 1 and refine mesh 3x (cells are 1/8 size).
            # Refining the output mesh is useful with a basis order of 2 or greater in the solution.
            output_basis_order = 1
            refine_levels = 3
        """
    }

    def __init__(self, name="outputsolndomain"):
        """Constructor.
        """
        OutputSoln.__init__(self, name)

    def preinitialize(self, problem):
        """Do mimimal initialization.
        """
        OutputSoln.preinitialize(self, problem)

        identifier = self.aliases[-1]
        self.writer.setFilename(problem.defaults.outputDir, problem.defaults.simName, identifier)

    def _configure(self):
        """Set members based using inventory.
        """
        OutputSoln._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleOutputSolnDomain.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """Factory associated with OutputSolnDomain.
    """
    return OutputSolnDomain()


# End of file
