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
# @file pythia.pyre/meshio/DataWriterHDF5.py
#
# @brief Python object for writing finite-element data to HDF5 file.

from .DataWriter import DataWriter
from .meshio import DataWriterHDF5 as ModuleDataWriterHDF5


class DataWriterHDF5(DataWriter, ModuleDataWriterHDF5):
    """Python object for writing finite-element data to HDF5 file.

    FACTORY: data_writer
    """

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="")
    filename.meta['tip'] = "Name of HDF5 file."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="datawriterhdf5"):
        """Constructor.
        """
        DataWriter.__init__(self, name)
        ModuleDataWriterHDF5.__init__(self)
        return

    def preinitialize(self):
        """Initialize writer.
        """
        DataWriter.preinitialize(self)
        return

    def setFilename(self, outputDir, simName, label):
        """Set filename from default options and inventory. If filename is given in inventory, use it,
        otherwise create filename from default options.
        """
        filename = self.filename or DataWriter.mkfilename(outputDir, simName, label, "h5")
        self.mkpath(filename)
        ModuleDataWriterHDF5.filename(self, filename)
        return

    # PRIVATE METHODS /////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleDataWriterHDF5.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////


def data_writer():
    """Factory associated with DataWriter.
    """
    return DataWriterHDF5()


# End of file
