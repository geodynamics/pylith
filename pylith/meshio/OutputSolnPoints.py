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
from .meshio import OutputSolnPoints as ModuleOutputSolnPoints


class OutputSolnPoints(OutputSoln, ModuleOutputSolnPoints):
    """
    Output of solution subfields at discrete points in the domain.

    :::{tip}
    Most output information can be configured at the problem level using the [`ProblemDefaults` Component](../problems/ProblemDefaults.md).
    :::

    Implements `OutputSoln`.
    """
    DOC_CONFIG = {
        "cfg": """
            [observer]
            label = stations
            data_fields = [displacement]

            # List of points where we want output.
            reader = pylith.meshio.PointsList
            reader.filename = stations.txt

            # Skip two time steps between output.
            output_trigger = pylith.meshio.OutputTriggerStep
            output_trigger.num_skip = 2

            # Write output to HDF5 file with name `domain.h5`.
            writer = pylith.meshio.DataWriterHDF5
            writer.filename = domain.h5

            output_basis_order = 1
        """
    }

    import pythia.pyre.inventory

    label = pythia.pyre.inventory.str("label", default="points")
    label.meta['tip'] = "Label identifier for points (used in constructing default filenames)."

    from .PointsList import PointsList
    reader = pythia.pyre.inventory.facility("reader", factory=PointsList, family="points_list")
    reader.meta['tip'] = "Reader for points list."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputsolnpoints"):
        """Constructor.
        """
        OutputSoln.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Do mimimal initialization.
        """
        OutputSoln.preinitialize(self, problem)

        stationNames, stationCoords = self.reader.read()

        # Convert to mesh coordinate system
        from spatialdata.geocoords.Converter import convert
        convert(stationCoords, problem.mesh().getCoordSys(), self.reader.coordsys)

        # Nondimensionalize
        stationCoords /= problem.normalizer.lengthScale.value

        ModuleOutputSolnPoints.setPoints(self, stationCoords, stationNames)

        identifier = self.aliases[-1]
        self.writer.setFilename(problem.defaults.outputDir, problem.defaults.simName, identifier)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleOutputSolnPoints.__init__(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def observer():
    """Factory associated with OutputSoln.
    """
    return OutputSolnPoints()


# End of file
