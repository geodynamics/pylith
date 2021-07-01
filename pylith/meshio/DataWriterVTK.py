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
# @file pythia.pyre/meshio/DataWriterVTK.py
#
# @brief Python object for writing finite-element data to VTK file.

from .DataWriter import DataWriter
from .meshio import DataWriterVTK as ModuleDataWriterVTK


class DataWriterVTK(DataWriter, ModuleDataWriterVTK):
    """Python object for writing finite-element data to VTK file.

    FACTORY: data_writer
    """

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="")
    filename.meta['tip'] = "Name of VTK file."

    timeFormat = pythia.pyre.inventory.str("time_format", default="%f")
    timeFormat.meta['tip'] = "C style format string for time stamp in filename."

    from pythia.pyre.units.time import second
    timeConstant = pythia.pyre.inventory.dimensional(
        "time_constant", default=1.0 * second, validator=pythia.pyre.inventory.greater(0.0 * second))
    timeConstant.meta['tip'] = "Values used to normalize time stamp in filename."

    precision = pythia.pyre.inventory.int("float_precision", default=6,
                                   validator=pythia.pyre.inventory.greater(0))
    precision.meta['tip'] = "Precision of floating point values in output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="datawritervtk"):
        """Constructor.
        """
        DataWriter.__init__(self, name)
        ModuleDataWriterVTK.__init__(self)
        return

    def preinitialize(self):
        """Initialize writer.
        """
        DataWriter.preinitialize(self)

        ModuleDataWriterVTK.timeFormat(self, self.timeFormat)
        ModuleDataWriterVTK.timeConstant(self, self.timeConstant.value)
        ModuleDataWriterVTK.precision(self, self.precision)
        return

    def setFilename(self, outputDir, simName, label):
        """Set filename from default options and inventory. If filename is given in inventory, use it,
        otherwise create filename from default options.
        """
        filename = self.filename or DataWriter.mkfilename(outputDir, simName, label, "vtk")
        self.mkpath(filename)
        ModuleDataWriterVTK.filename(self, filename)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Configure object.
        """
        DataWriter._configure(self)
        return

    # PRIVATE METHODS /////////////////////////////////////////////////////

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
