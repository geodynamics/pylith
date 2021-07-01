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
# @file pythia.pyre/meshio/OutputSolnPoints.py
#
# @brief Python object for managing output of finite-element solution
# information over a subdomain.
#
# FACTORY: observer

from .OutputSoln import OutputSoln
from .meshio import OutputSolnPoints as ModuleOutputSolnPoints


class OutputSolnPoints(OutputSoln, ModuleOutputSolnPoints):
    """Python object for managing output of finite-element solution
    information over a subdomain.

    FACTORY: observer
    """

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
