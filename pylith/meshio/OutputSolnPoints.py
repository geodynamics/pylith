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
# @file pyre/meshio/OutputSolnPoints.py
#
# @brief Python object for managing output of finite-element solution
# information over a subdomain.
#
# FACTORY: observer

from OutputManager import OutputManager
from meshio import OutputSolnPoints as ModuleOutputSolnPoints


def validateFilename(value):
    """
    Validate filename with list of points.
    """
    if 0 == len(value):
        raise ValueError("Filename for list of points not specified.")
    return value


class OutputSolnPoints(OutputManager, ModuleOutputSolnPoints):
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

    import pyre.inventory

    from PointsList import PointsList
    reader = pyre.inventory.facility("reader", factory=PointsList, family="points_list")
    reader.meta['tip'] = "Reader for points list."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputsolnpoints"):
        """
        Constructor.
        """
        OutputManager.__init__(self, name)
        return

    def preinitialize(self):
        """
        Do
        """
        OutputManager.preinitialize(self)
        return

    def initialize(self, mesh, normalizer):
        """
        Initialize output manager.
        """
        logEvent = "%sinit" % self._loggingPrefix
        self._eventLogger.eventBegin(logEvent)

        OutputManager.initialize(self, normalizer)

        # Read points
        stations, points = self.reader.read()

        # Convert to mesh coordinate system
        from spatialdata.geocoords.Converter import convert
        convert(points, mesh.coordsys(), self.coordsys)

        ModuleOutputSolnPoints.setupInterpolator(self, mesh, points, stations, normalizer)
        self.mesh = ModuleOutputSolnPoints.pointsMesh(self)

        self._eventLogger.eventEnd(logEvent)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        OutputManager._configure(self)
        return

    def _createModuleObj(self):
        """
        Create handle to C++ object.
        """
        ModuleOutputSolnPoints.__init__(self)
        return

    def _open(self, mesh, nsteps, label, labelId):
        """
        Call C++ open();
        """
        if label != None and labelId != None:
            ModuleOutputSolnPoints.open(self, mesh, nsteps, label, labelId)
        else:
            ModuleOutputSolnPoints.open(self, mesh, nsteps)

        ModuleOutputSolnPoints.writePointNames(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """
    Factory associated with OutputManager.
    """
    return OutputSolnPoints()


# End of file
