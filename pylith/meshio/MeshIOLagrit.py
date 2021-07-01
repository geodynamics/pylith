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
# @file pythia.pyre/meshio/MeshIOLagrit.py
#
# @brief Python object for reading/writing finite-element mesh from
# LaGriT.
#
# Factory: mesh_io

from .MeshIOObj import MeshIOObj
from .meshio import MeshIOLagrit as ModuleMeshIOLagrit


def validateFilenameGmv(value):
    """Validate filename.
    """
    if 0 == len(value):
        raise ValueError("Filename for LaGriT input mesh not specified.")
    try:
        open(value, "r")
    except IOError:
        raise IOError("LaGriT input mesh '{}' not found.".format(value))
    return value


def validateFilenamePset(value):
    """Validate filename.
    """
    if 0 == len(value):
        raise ValueError("Filename for LaGriT pset file not specified.")
    try:
        open(value, "r")
    except IOError:
        raise IOError("LaGriT pset file '{}' not found.".format(value))
    return value


class MeshIOLagrit(MeshIOObj, ModuleMeshIOLagrit):
    """Python object for reading/writing finite-element mesh from LaGriT.

    FACTORY: mesh_io
    """

    import pythia.pyre.inventory

    filenameGmv = pythia.pyre.inventory.str("filename_gmv", default="mesh.gmv",
                                     validator=validateFilenameGmv)
    filenameGmv.meta['tip'] = "Name of mesh GMV file."

    filenamePset = pythia.pyre.inventory.str("filename_pset", default="mesh.pset",
                                      validator=validateFilenamePset)
    filenamePset.meta['tip'] = "Name of mesh PSET file."

    flipEndian = pythia.pyre.inventory.bool("flip_endian", default=False)
    flipEndian.meta['tip'] = "Flip endian type when reading/writing binary " \
                             "files."

    ioInt32 = pythia.pyre.inventory.bool("io_int32", default=True)
    ioInt32.meta['tip'] = "PSTE files use 32-bit integers."

    isRecordHeader32Bit = pythia.pyre.inventory.bool("record_header_32bit",
                                              default=True)
    isRecordHeader32Bit.meta['tip'] = "Fortran record header is 32-bit."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys",
                                       factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="meshiolagrit"):
        """Constructor.
        """
        MeshIOObj.__init__(self, name)
        return

    def preinitialize(self):
        """Do minimal initialization."""
        MeshIOObj.preinitialize(self)

        ModuleMeshIOLagrit.filenameGmv(self, self.filenameGmv)
        ModuleMeshIOLagrit.filenamePset(self, self.filenamePset)
        ModuleMeshIOLagrit.flipEndian(self, self.flipEndian)
        ModuleMeshIOLagrit.ioInt32(self, self.ioInt32)
        ModuleMeshIOLagrit.isRecordHeader32Bit(self, self.isRecordHeader32Bit)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        MeshIOObj._configure(self)
        return

    def _createModuleObj(self):
        """Create C++ MeshIOLagrit object.
        """
        ModuleMeshIOLagrit.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
    """Factory associated with MeshIOLagrit.
    """
    return MeshIOLagrit()


# End of file
