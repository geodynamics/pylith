# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
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
    """
    Python object for managing output of finite-element solution
    information over a subdomain.

    INVENTORY

    Properties
      - None

    Facilities
      - *reader* Reader for list of points.

    FACTORY: observer
    """

    # INVENTORY //////////////////////////////////////////////////////////

    import pythia.pyre.inventory

    label = pythia.pyre.inventory.str("label", default="points")
    label.meta['tip'] = "Label identifier for points (used in constructing default filenames)."

    from .PointsList import PointsList
    reader = pythia.pyre.inventory.facility("reader", factory=PointsList, family="points_list")
    reader.meta['tip'] = "Reader for points list."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputsolnpoints"):
        """
        Constructor.
        """
        OutputSoln.__init__(self, name)
        return

    def preinitialize(self, problem):
        """
        Do mimimal initialization.
        """
        OutputSoln.preinitialize(self, problem)

        stationNames, stationCoords = self.reader.read()

        # Convert to mesh coordinate system
        from spatialdata.geocoords.Converter import convert
        convert(points, problem.mesh.coordsys(), reader.coordsys)

        # Nondimensionalize
        stationsCoords /= problem.normalizer.lengthScale.value

        ModuleOutputSolnPoints.stations(stationCoords, stationNames)

        identifier = self.aliases[-1]
        self.writer.setFilename(problem.defaults.outputDir, problem.defaults.simName, identifier)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        OutputSoln._configure(self)
        return

    def _createModuleObj(self):
        """
        Create handle to C++ object.
        """
        ModuleOutputSolnPoints.__init__(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def observer():
    """
    Factory associated with OutputSoln.
    """
    return OutputSolnPoints()


# End of file
