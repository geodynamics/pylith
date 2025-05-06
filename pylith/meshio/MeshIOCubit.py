# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file pythia.pyre/meshio/MeshIOCubit.py
#
# @brief Python object for reading/writing finite-element mesh from
# Cubit.
#
# Factory: mesh_input, mesh_output

import pathlib

from .MeshIOObj import MeshIOObj
from .meshio import MeshIOCubit as ModuleMeshIOCubit


class MeshIOCubit(MeshIOObj, ModuleMeshIOCubit):
    """
    Reader for finite-element meshes from Exodus II files (usually from Cubit).

    :::{warning}
    The coordinate system associated with the mesh must be a Cartesian coordinate system, such as a generic Cartesian coordinate system or a geographic projection.
    :::

    Implements `MeshIOObj`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.mesh_generator.reader]
            filename = mesh_quad.exo
            coordsys.space_dim = 2
        """
    }

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="mesh.exo")
    filename.meta['tip'] = "Name of Cubit Exodus file."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys", factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."

    def __init__(self, name="meshiocubit"):
        """Constructor.
        """
        mode = MeshIOObj.READ
        MeshIOObj.__init__(self, mode, name)

    def preinitialize(self):
        """Do minimal initialization."""
        MeshIOObj.preinitialize(self)
        ModuleMeshIOCubit.setFilename(self, self.filename)

    def _configure(self):
        """Set members based using inventory.
        """
        MeshIOObj._configure(self)

    def _createModuleObj(self):
        """Create C++ MeshIOCubit object.
        """
        ModuleMeshIOCubit.__init__(self)

    def _validate(self, context):
        if 0 == len(self.filename):
            context.error(ValueError("Filename for CUBIT mesh not specified."))
        if self.mode == self.READ and not pathlib.Path(self.filename).is_file():
            context.error(IOError(f"CUBIT input mesh '{self.filename}' not found."))


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_input():
    """Factory associated with MeshIOCubit.
    """
    return MeshIOCubit()


# End of file
