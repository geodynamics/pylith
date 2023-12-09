# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

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

    labelName = pythia.pyre.inventory.str("label", default="", validator=validateLabel)
    labelName.meta['tip'] = "Name of label identifier for external boundary."

    labelValue = pythia.pyre.inventory.int("label_value", default=1)
    labelValue.meta['tip'] = "Value of label identifier for external boundary (tag of physical group in Gmsh files)."

    def __init__(self, name="outputsolnsubset"):
        """Constructor.
        """
        OutputSoln.__init__(self, name)

    def preinitialize(self, problem):
        """Do mimimal initialization.
        """
        OutputSoln.preinitialize(self, problem)
        ModuleOutputSolnBoundary.setLabelName(self, self.labelName)
        ModuleOutputSolnBoundary.setLabelValue(self, self.labelValue)

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
