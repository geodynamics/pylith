# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

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
    """
    Reader for finite-element meshes using a simple ASCII format.

    :::{warning}
    The coordinate system associated with the mesh must be a Cartesian coordinate system, such as a generic Cartesian coordinate system or a geographic projection.
    :::

    Implements `MeshIOObj`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.mesh_generator.reader]
            filename = mesh_quad.txt
            coordsys.space_dim = 2
        """
    }

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="", validator=validateFilename)
    filename.meta['tip'] = "Name of mesh file"

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys", factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."

    def __init__(self, name="meshioascii"):
        """Constructor.
        """
        MeshIOObj.__init__(self, name)

    def preinitialize(self):
        """Do minimal initialization."""
        MeshIOObj.preinitialize(self)
        ModuleMeshIOAscii.setFilename(self, self.filename)

    def _configure(self):
        """Set members based using inventory.
        """
        MeshIOObj._configure(self)

    def _createModuleObj(self):
        """Create C++ MeshIOAscii object.
        """
        ModuleMeshIOAscii.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
    """Factory associated with MeshIOAscii.
    """
    return MeshIOAscii()


# End of file
