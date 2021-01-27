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
# @file pyre/meshio/DataWriterVTK.py
#
# @brief Python object for writing finite-element data to VTK file.

from .DataWriter import DataWriter
from .meshio import DataWriterVTK as ModuleDataWriterVTK


class DataWriterVTK(DataWriter, ModuleDataWriterVTK):
    """
    Python object for writing finite-element data to VTK file.

    INVENTORY

    Properties
      - *filename* Name of VTK file.
      - *time_format* C style format string for time stamp in filename.
      - *time_constant* Value used to normalize time stamp in filename.

    Facilities
      - None

    FACTORY: data_writer
    """

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="")
    filename.meta['tip'] = "Name of VTK file."

    timeFormat = pyre.inventory.str("time_format", default="%f")
    timeFormat.meta['tip'] = "C style format string for time stamp in filename."

    from pyre.units.time import second
    timeConstant = pyre.inventory.dimensional(
        "time_constant", default=1.0 * second, validator=pyre.inventory.greater(0.0 * second))
    timeConstant.meta['tip'] = "Values used to normalize time stamp in filename."

    precision = pyre.inventory.int("float_precision", default=6,
                                   validator=pyre.inventory.greater(0))
    precision.meta['tip'] = "Precision of floating point values in output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="datawritervtk"):
        """
        Constructor.
        """
        DataWriter.__init__(self, name)
        ModuleDataWriterVTK.__init__(self)
        return

    def preinitialize(self):
        """
        Initialize writer.
        """
        DataWriter.preinitialize(self)

        ModuleDataWriterVTK.timeFormat(self, self.timeFormat)
        ModuleDataWriterVTK.timeConstant(self, self.timeConstant.value)
        ModuleDataWriterVTK.precision(self, self.precision)
        return

    def setFilename(self, outputDir, simName, label):
        """
        Set filename from default options and inventory. If filename is given in inventory, use it,
        otherwise create filename from default options.
        """
        filename = self.filename or DataWriter.mkfilename(outputDir, simName, label, "vtk")
        self.mkpath(filename)
        ModuleDataWriterVTK.filename(self, filename)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Configure object.
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
    """
    Factory associated with DataWriter.
    """
    return DataWriterVTK()


# End of file
