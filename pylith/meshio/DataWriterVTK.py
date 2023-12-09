# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .DataWriter import DataWriter
from .meshio import DataWriterVTK as ModuleDataWriterVTK


class DataWriterVTK(DataWriter, ModuleDataWriterVTK):
    """
    Writer of solution, auxiliary, and derived subfields to a VTK file.

    Implements `DataWriter`.
    """
    DOC_CONFIG = {
        "cfg": """
            [data_writer]
            filename = domain_solution.vtk
            time_format = %0.2f
            time_constant = 1.0*year
            float_precision = 6
        """
    }

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="")
    filename.meta['tip'] = "Name of VTK file."

    timeFormat = pythia.pyre.inventory.str("time_format", default="%f")
    timeFormat.meta['tip'] = "C style format string for time stamp in filename."

    from pythia.pyre.units.time import second
    timeConstant = pythia.pyre.inventory.dimensional("time_constant", default=1.0 * second, validator=pythia.pyre.inventory.greater(0.0 * second))
    timeConstant.meta['tip'] = "Values used to normalize time stamp in filename."

    precision = pythia.pyre.inventory.int("float_precision", default=6, validator=pythia.pyre.inventory.greater(0))
    precision.meta['tip'] = "Precision of floating point values in output."

    def __init__(self, name="datawritervtk"):
        """Constructor.
        """
        DataWriter.__init__(self, name)
        ModuleDataWriterVTK.__init__(self)

    def preinitialize(self):
        """Initialize writer.
        """
        DataWriter.preinitialize(self)

        ModuleDataWriterVTK.timeFormat(self, self.timeFormat)
        ModuleDataWriterVTK.timeConstant(self, self.timeConstant.value)
        ModuleDataWriterVTK.precision(self, self.precision)

    def setFilename(self, outputDir, simName, label):
        """Set filename from default options and inventory. If filename is given in inventory, use it,
        otherwise create filename from default options.
        """
        filename = self.filename or DataWriter.mkfilename(outputDir, simName, label, "vtk")
        self.mkpath(filename)
        ModuleDataWriterVTK.filename(self, filename)

    def _configure(self):
        """Configure object.
        """
        DataWriter._configure(self)

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleDataWriterVTK.__init__(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def data_writer():
    """Factory associated with DataWriter.
    """
    return DataWriterVTK()


# End of file
