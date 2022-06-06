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
from .meshio import MeshIOPetsc as ModuleMeshIOPetsc


class MeshIOPetsc(MeshIOObj, ModuleMeshIOPetsc):
    """
    Python object for reading/writing finite-element mesh from simple PETSc file.
    """

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="")
    filename.meta['tip'] = "Name of mesh file for reading with PETSc."

    prefix = pythia.pyre.inventory.str("options_prefix", default="")
    prefix.meta['tip'] = "Name of PETSc options prefix for this mesh."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys", factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."

    def __init__(self, name="meshiopetsc"):
        """Constructor.
        """
        MeshIOObj.__init__(self, name)

    def preinitialize(self):
        """Do minimal initialization."""
        MeshIOObj.preinitialize(self)
        ModuleMeshIOPetsc.setFilename(self, self.filename)
        ModuleMeshIOPetsc.setPrefix(self, self.prefix)

    def _configure(self):
        """Set members based using inventory.
        """
        MeshIOObj._configure(self)

    def _createModuleObj(self):
        """Create C++ MeshIOPetsc object.
        """
        ModuleMeshIOPetsc.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
    """Factory associated with MeshIOPetsc.
    """
    return MeshIOPetsc()


# End of file
