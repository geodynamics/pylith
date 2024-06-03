# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file pythia.pyre/meshio/MeshIOCubit.py
#
# @brief Python object for reading/writing finite-element mesh from
# Cubit.
#
# Factory: mesh_io

from .MeshIOObj import MeshIOObj
from .meshio import MeshIOCubit as ModuleMeshIOCubit


def validateFilename(value):
    """Validate filename.
    """
    if 0 == len(value):
        raise ValueError("Filename for CUBIT/Trlis input mesh not specified.")
    try:
        open(value, "r")
    except IOError:
        raise IOError("CUBIT/Trelis input mesh '{}' not found.".format(value))
    return value


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
            use_nodeset_names = True
            coordsys.space_dim = 2
        """
    }

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="mesh.exo", validator=validateFilename)
    filename.meta['tip'] = "Name of Cubit Exodus file."

    useNames = pythia.pyre.inventory.bool("use_nodeset_names", default=True)
    useNames.meta['tip'] = "Use nodeset names instead of ids."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys", factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."

    def __init__(self, name="meshiocubit"):
        """Constructor.
        """
        MeshIOObj.__init__(self, name)

    def preinitialize(self):
        """Do minimal initialization."""
        MeshIOObj.preinitialize(self)
        ModuleMeshIOCubit.setFilename(self, self.filename)
        ModuleMeshIOCubit.setUseNodesetNames(self, self.useNames)

    def _configure(self):
        """Set members based using inventory.
        """
        MeshIOObj._configure(self)

    def _createModuleObj(self):
        """Create C++ MeshIOCubit object.
        """
        ModuleMeshIOCubit.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
    """Factory associated with MeshIOCubit.
    """
    return MeshIOCubit()


# End of file
