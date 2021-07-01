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
# @file pythia.pyre/meshio/MeshIOAscii.py
#
# @brief Python object for reading/writing finite-element mesh from
# simple ASCII file.
#
# Factory: mesh_io

from .MeshIOObj import MeshIOObj
from .meshio import MeshIOAscii as ModuleMeshIOAscii


def validateFilename(value):
    """Validate filename.
    """
    if 0 == len(value):
        msg = "Filename for ASCII input mesh not specified.  " + \
            "To test PyLith, run an example as discussed in the manual."
        raise ValueError(msg)
    try:
        open(value, "r")
    except IOError:
        raise IOError("ASCII input mesh '{}' not found.".format(value))
    return value


class MeshIOAscii(MeshIOObj, ModuleMeshIOAscii):
    """Python object for reading/writing finite-element mesh from simple
    ASCII file.

    Factory: mesh_io
    """

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="",
                                  validator=validateFilename)
    filename.meta['tip'] = "Name of mesh file"

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys",
                                       factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="meshioascii"):
        """Constructor.
        """
        MeshIOObj.__init__(self, name)
        return

    def preinitialize(self):
        """Do minimal initialization."""
        MeshIOObj.preinitialize(self)

        ModuleMeshIOAscii.filename(self, self.filename)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        MeshIOObj._configure(self)
        return

    def _createModuleObj(self):
        """Create C++ MeshIOAscii object.
        """
        ModuleMeshIOAscii.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
    """Factory associated with MeshIOAscii.
    """
    return MeshIOAscii()


# End of file
