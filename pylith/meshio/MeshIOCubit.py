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
    """Python object for reading/writing finite-element mesh from Cubit.

    FACTORY: mesh_io
    """

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="mesh.exo", validator=validateFilename)
    filename.meta['tip'] = "Name of Cubit Exodus file."

    useNames = pythia.pyre.inventory.bool("use_nodeset_names", default=True)
    useNames.meta['tip'] = "Use nodeset names instead of ids."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys",
                                       factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="meshiocubit"):
        """Constructor.
        """
        MeshIOObj.__init__(self, name)
        return

    def preinitialize(self):
        """Do minimal initialization."""
        MeshIOObj.preinitialize(self)

        ModuleMeshIOCubit.filename(self, self.inventory.filename)
        ModuleMeshIOCubit.useNodesetNames(self, self.inventory.useNames)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        MeshIOObj._configure(self)
        return

    def _createModuleObj(self):
        """Create C++ MeshIOCubit object.
        """
        ModuleMeshIOCubit.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
    """Factory associated with MeshIOCubit.
    """
    return MeshIOCubit()


# End of file
