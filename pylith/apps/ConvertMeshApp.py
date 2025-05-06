# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# Application for converting mesh files from one format to another.

import pathlib

from .PetscApplication import PetscApplication

def validateFilename(value):
    """Validate filename.
    """
    if 0 == len(value):
        msg = "Filename for input mesh not specified."
        raise ValueError(msg)
    if not pathlib.Path(value).is_file():
        raise IOError(f"Input mesh '{value}' not found.")
    return value


class ConvertMeshApp(PetscApplication):
    """Application for converting mesh files from one format to another.
    """

    import pythia.pyre.inventory
    from pylith.meshio.MeshIOPetsc import MeshIOPetsc

    reader = pythia.pyre.inventory.facility("reader", family="mesh_input", factory=MeshIOPetsc)
    reader.meta["tip"] = "Reader for input mesh file."

    writer = pythia.pyre.inventory.facility("writer", family="mesh_output", factory=MeshIOPetsc)
    writer.meta["tip"] = "Writer for output mesh file."

    checkTopology = pythia.pyre.inventory.bool("check_topology", default=True)
    checkTopology.meta['tip'] = "Check topology of imported mesh."

    def __init__(self, name="convertmeshapp"):
        """Constructor.
        """
        PetscApplication.__init__(self, name)

    def main(self, *args, **kwds):
        """Run the application.
        """
        from pylith.meshio.meshio import MeshConverter

        self.initialize()
        MeshConverter.convert(self.writer, self.reader, self.checkTopology)
                
    def _configure(self):
        """Setup members using inventory.
        """
        PetscApplication._configure(self)

    def initialize(self):
        self.reader.preinitialize()
        self.writer.preinitialize()

# End of file
