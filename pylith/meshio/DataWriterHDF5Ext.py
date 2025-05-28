# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .DataWriter import DataWriter
from .meshio import DataWriterHDF5Ext as ModuleDataWriterHDF5Ext


class DataWriterHDF5Ext(DataWriter, ModuleDataWriterHDF5Ext):
    """
    Writer of solution, auxiliary, and derived subfields to an HDF5 file with datasets stored in external binary files.

    Implements `DataWriter`.
    """
    DOC_CONFIG = {
        "cfg": """
            [data_writer]
            filename = domain_solution.h5
        """
    }

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="")
    filename.meta['tip'] = "Name of HDF5 file."

    def __init__(self, name="datawriterhdf5"):
        """Constructor.
        """
        DataWriter.__init__(self, name)

    def preinitialize(self):
        """Initialize writer.
        """
        DataWriter.preinitialize(self)

    def setFilename(self, outputDir, simName, label):
        """Set filename from default options and inventory. If filename is given in inventory, use it,
        otherwise create filename from default options.
        """
        filename = self.filename or DataWriter.makeFilename(outputDir, simName, label, "h5")
        self.makePath(filename)
        ModuleDataWriterHDF5Ext.filename(self, str(filename))

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

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleDataWriterHDF5Ext.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def data_writer():
    """Factory associated with DataWriter.
    """
    return DataWriterHDF5Ext()


# End of file
