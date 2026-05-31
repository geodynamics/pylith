# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pythia.pyre.components.Component import Component

from .Initializer import Initializer as InitializerBase
from .Initializer import phaseFactory


class Convert(Component):
    """
    Mesh initialization phases for reordering and converting to a different format.

    :::{seealso}
    See [`Initializer` Component](Initializer.md).
    :::
    """

    import pythia.pyre.inventory

    from .MeshReader import MeshReader
    from .MeshReordering import MeshReordering
    from .MeshWriter import MeshWriter

    read_mesh = pythia.pyre.inventory.facility(
        "read_mesh", family="initialize_phase", factory=MeshReader
    )
    read_mesh.meta["tip"] = "Read mesh."

    reorder_mesh = pythia.pyre.inventory.facility(
        "reorder_mesh", family="initialize_phase", factory=MeshReordering
    )
    reorder_mesh.meta["tip"] = "Reorder mesh using reverse Cuthill-McKee algorithm."

    write_mesh = pythia.pyre.inventory.facility(
        "write_mesh", family="initialize_phase", factory=MeshWriter
    )
    write_mesh.meta["tip"] = "Write mesh."


    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="convert_phases"):
        """Constructor."""
        Component.__init__(self, name, facility="convert_phases")

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to ensure order is [displacement, pressure, trace_strain].

        """
        return [
            self.read_mesh,
            self.reorder_mesh,
            self.write_mesh,
        ]



class Initializer(InitializerBase):
    """Python mesh initializer with convert phases.
    """

    import pythia.pyre.inventory

    phases = pythia.pyre.inventory.facilityArray("phases", family="convert_phases", itemFactory=phaseFactory, factory=Convert)
    phases.meta['tip'] = "Phases in converting mesh."


# FACTORIES ////////////////////////////////////////////////////////////
def mesh_initializer():
    """Factory associated with Initializer.
    """
    return Initializer()


# End of file
