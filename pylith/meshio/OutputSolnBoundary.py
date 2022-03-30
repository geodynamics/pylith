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

from .OutputSoln import OutputSoln
from .meshio import OutputSolnBoundary as ModuleOutputSolnBoundary


def validateLabel(value):
    """Validate label for group/nodeset/pset.
    """
    if not value.strip():
        raise ValueError("Label for group/nodeset/pset in mesh not specified.")
    return value


class OutputSolnBoundary(OutputSoln, ModuleOutputSolnBoundary):
    """
    Output of solution subfields over an external boundary.

    :::{tip}
    Most output information can be configured at the problem level using the [`ProblemDefaults` Component](../problems/ProblemDefaults.md).
    :::

    Implements `OutputSoln`.
    """
    DOC_CONFIG = {
        "cfg": """
            [observer]
            data_fields = [displacement]

            label = boundary_xpos

            # Skip two time steps between output.
            output_trigger = pylith.meshio.OutputTriggerStep
            output_trigger.num_skip = 2

            # Write output to HDF5 file with name `boundary_xpos.h5`.
            writer = pylith.meshio.DataWriterHDF5
            writer.filename = boundary_xpos.h5

            output_basis_order = 1
        """
    }

    import pythia.pyre.inventory

    label = pythia.pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for external boundary."

    def __init__(self, name="outputsolnsubset"):
        """Constructor.
        """
        OutputSoln.__init__(self, name)

    def preinitialize(self, problem):
        """Do mimimal initialization.
        """
        OutputSoln.preinitialize(self, problem)
        ModuleOutputSolnBoundary.setLabel(self, self.label)

        identifier = self.aliases[-1]
        self.writer.setFilename(problem.defaults.outputDir, problem.defaults.simName, identifier)

    def _configure(self):
        """Set members based using inventory.
        """
        OutputSoln._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleOutputSolnBoundary.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """Factory associated with OutputSoln.
    """
    return OutputSolnBoundary()


# End of file
