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
# @file pythia.pyre/meshio/DataWriterHDF5Ext.py
#
# @brief Python object for writing finite-element data to HDF5 file
# with datasets stored in external binary files.

from .DataWriter import DataWriter
from .meshio import DataWriterHDF5Ext as ModuleDataWriterHDF5Ext


class DataWriterHDF5Ext(DataWriter, ModuleDataWriterHDF5Ext):
    """
    @brief Python object for writing finite-element data to HDF5 file
    with datasets stored in external binary files.

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
        ModuleDataWriterHDF5Ext.filename(self, filename)
        return

    def close(self):
        """Close writer.
        """
        ModuleDataWriterHDF5Ext.close(self)

        # Only write Xdmf file on proc 0
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if not comm.rank:
            from .Xdmf import Xdmf
            xdmf = Xdmf()
            xdmf.write(ModuleDataWriterHDF5Ext.hdf5Filename(
                self), verbose=False)
        return

    # PRIVATE METHODS /////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleDataWriterHDF5Ext.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def data_writer():
    """Factory associated with DataWriter.
    """
    return DataWriterHDF5Ext()


# End of file
