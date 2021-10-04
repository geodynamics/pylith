# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pythia.pyre/meshio/MeshIOPETSc.py
#
# @brief Python object for reading/writing finite-element mesh from
# simple gmsh file.
#
# Factory: mesh_io

from .MeshIOObj import MeshIOObj
from .meshio import MeshIOPETSc as ModuleMeshIOPETSc


class MeshIOPETSc(MeshIOObj, ModuleMeshIOPETSc):
    """Python object for reading/writing finite-element mesh from simple
    PETSc file.

    Factory: mesh_io
    """

    import pythia.pyre.inventory

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys",
                                       factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="meshiopetsc"):
        """Constructor.
        """
        MeshIOObj.__init__(self, name)
        return

    def preinitialize(self):
        """Do minimal initialization."""
        MeshIOObj.preinitialize(self)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        MeshIOObj._configure(self)
        return

    def _createModuleObj(self):
        """Create C++ MeshIOPETSc object.
        """
        ModuleMeshIOPETSc.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
    """Factory associated with MeshIOPETSc.
    """
    return MeshIOPETSc()


# End of file
