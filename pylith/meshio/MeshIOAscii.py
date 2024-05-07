# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import pathlib

from .MeshIOObj import MeshIOObj
from .meshio import MeshIOAscii as ModuleMeshIOAscii


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

    filename = pythia.pyre.inventory.str("filename", default="")
    filename.meta['tip'] = "Name of mesh file"

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys", factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."

    def __init__(self, mode=MeshIOObj.READ, name="meshioascii"):
        """Constructor.
        """
        MeshIOObj.__init__(self, mode, name)

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

    def _validate(self, context):
        if 0 == len(self.filename):
            context.error(ValueError("Filename for ASCII mesh not specified."))
        if self.mode == self.READ and not pathlib.Path(self.filename).is_file():
            context.error(IOError(f"ASCII input mesh '{self.filename}' not found."))

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_input():
    """Factory associated with MeshIOAscii.
    """
    return MeshIOAscii(mode=MeshIOAscii.READ)


def mesh_output():
    """Factory associated with MeshIOAscii.
    """
    return MeshIOAscii(mode=MeshIOAscii.WRITE)


# End of file
