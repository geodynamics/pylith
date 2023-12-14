# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file pythia.pyre/meshio/OutputPhysics.py
#
# @brief Python object for managing output over points with constrained degrees of freedom.
#
# Factory: output_manager

from .OutputObserver import OutputObserver
from .meshio import OutputPhysics as ModuleOutputPhysics


class OutputPhysics(OutputObserver, ModuleOutputPhysics):
    """
    Output for objects implementing physics (materials and boundary conditions).

    :::{tip}
    Most output information can be configured at the problem level using the [`ProblemDefaults` Component](../problems/ProblemDefaults.md).
    :::

    Implements `OutputObserver`.
    """
    DOC_CONFIG = {
        "cfg": """
            [observer]
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

    infoFields = pythia.pyre.inventory.list("info_fields", default=["all"])
    infoFields.meta['tip'] = "Names of auxiliary subfields to include in info output."

    dataFields = pythia.pyre.inventory.list("data_fields", default=["all"])
    dataFields.meta['tip'] = "Names of solution, auxiliary, and derived subfields to include in data output."

    def __init__(self, name="outputphysics"):
        """Constructor.
        """
        OutputObserver.__init__(self, name)

    def preinitialize(self, problem, identifier):
        """Do mimimal initialization.
        """
        OutputObserver.preinitialize(self, problem)
        ModuleOutputPhysics.setInfoFields(self, self.infoFields)
        ModuleOutputPhysics.setDataFields(self, self.dataFields)

        self.writer.setFilename(problem.defaults.outputDir, problem.defaults.simName, identifier)

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleOutputPhysics.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """Factory associated with OutputObserver.
    """
    return OutputPhysics()


# End of file
