# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/sources/Source.py
#
# @brief Python abstract base class for managing input and out put 
# sources not necessarily pertaining to domain boundaries
#
# Factory: source

from pylith.problems.Physics import Physics
from .sources import Source as ModuleSource


def validateLabel(value):
    """Validate descriptive label.
    """
    if 0 == len(value):
        raise ValueError("Descriptive label for source not specified.")
    return value


class Source(Physics, ModuleSource):
    """Python source property manager.

    FACTORY: source
    """

    import pythia.pyre.inventory

    sourceId = pythia.pyre.inventory.int("id", default=0)
    sourceId.meta['tip'] = "Source identifier (from mesh generator)."

    label = pythia.pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Descriptive label for source."

    from pylith.meshio.PointsList import PointsList
    reader = pythia.pyre.inventory.facility("reader", factory=PointsList, family="points_list")
    reader.meta['tip'] = "Reader for points list."

    def __init__(self, name="source"):
        """Constructor.
        """
        Physics.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Setup source.
        """
        Physics.preinitialize(self, problem)

        ModuleSource.setSourceId(self, self.sourceId)
        ModuleSource.setDescriptiveLabel(self, self.label)

        sourceNames, sourceCoords = self.reader.read()

        # Convert to mesh coordinate system
        from spatialdata.geocoords.Converter import convert
        convert(sourceCoords, problem.mesh().getCoordSys(), self.reader.coordsys)

        # Nondimensionalize
        sourceCoords /= problem.normalizer.lengthScale.value

        ModuleSource.setPoints(self, sourceCoords, sourceNames)
        return


# End of file
