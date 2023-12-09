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
from .meshio import DataWriterHDF5 as ModuleDataWriterHDF5


class DataWriterHDF5(DataWriter, ModuleDataWriterHDF5):
    """
    Writer of solution, auxiliary, and derived subfields to an HDF5 file.

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
        ModuleDataWriterHDF5.__init__(self)

    def preinitialize(self):
        """Initialize writer.
        """
        DataWriter.preinitialize(self)

    def setFilename(self, outputDir, simName, label):
        """Set filename from default options and inventory. If filename is given in inventory, use it,
        otherwise create filename from default options.
        """
        filename = self.filename or DataWriter.mkfilename(outputDir, simName, label, "h5")
        self.mkpath(filename)
        ModuleDataWriterHDF5.filename(self, filename)

    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleDataWriterHDF5.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////


def data_writer():
    """Factory associated with DataWriter.
    """
    return DataWriterHDF5()


# End of file
