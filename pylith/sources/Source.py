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


def validateDescription(value):
    """Validate description.
    """
    if 0 == len(value):
        raise ValueError("Description for material not specified.")
    return value


class Source(Physics, ModuleSource):
    """Python source property manager.

    FACTORY: source
    """

    import pythia.pyre.inventory

    description = pythia.pyre.inventory.str("description", default="", validator=validateDescription)
    description.meta['tip'] = "Descriptive label for material."

    labelName = pythia.pyre.inventory.str("label", default="source-id", validator=pythia.pyre.inventory.choice(["source-id"]))
    labelName.meta['tip'] = "Name of label for source. Currently only 'source-id' is allowed."

    labelValue = pythia.pyre.inventory.int("label_value", default=1)
    labelValue.meta['tip'] = "Value of label identifying source."

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
        ModuleSource.setDescription(self, self.description)
        ModuleSource.setLabelName(self, self.labelName)
        ModuleSource.setLabelValue(self, self.labelValue)

        sourceNames, sourceCoords = self.reader.read()

        # Convert to mesh coordinate system
        from spatialdata.geocoords.Converter import convert
        convert(sourceCoords, problem.mesh().getCoordSys(), self.reader.coordsys)

        # Nondimensionalize
        if hasattr(problem.normalizer,'lengthScale'):
            sourceCoords /= problem.normalizer.lengthScale.value
        else:
            sourceCoords /= (problem.normalizer.shearWaveSpeed.value * problem.normalizer.wavePeriod.value)

        ModuleSource.setPoints(self, sourceCoords, sourceNames)
        return


# End of file
