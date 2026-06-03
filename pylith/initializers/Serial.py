# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pythia.pyre.components.Component import Component


class Serial(Component):
    """
    Mesh initialization phases for reading in serial.

    :::{seealso}
    See [`Initializer` Component](Initializer.md).
    :::
    """

    DOC_CONFIG = {
        "cfg": """
            # Equivalent manual construction
            phases = [read_mesh, distribute_mesh, insert_interfaces, refine_mesh]
            read_mesh = pylith.initializers.MeshReader
            distribute_mesh = pylith.initializers.MeshDistributor
            insert_interfaces = pylith.initializers.InsertInterfaces
            refine_mesh = pylith.initializers.MeshRefiner
        """
    }

    import pythia.pyre.inventory

    from .MeshReader import MeshReader
    from .MeshReordering import MeshReordering
    from .MeshInsertInterfaces import MeshInsertInterfaces
    from .MeshDistributor import MeshDistributor
    from .MeshRefiner import MeshRefiner

    read_mesh = pythia.pyre.inventory.facility(
        "read_mesh", family="initialize_phase", factory=MeshReader
    )
    read_mesh.meta["tip"] = "Read mesh in serial."

    reorder_mesh = pythia.pyre.inventory.facility(
        "reorder_mesh", family="initialize_phase", factory=MeshReordering
    )
    reorder_mesh.meta["tip"] = "Reorder mesh using reverse Cuthill-McKee algorithm."

    distribute_mesh = pythia.pyre.inventory.facility(
        "distribute_mesh", family="initialize_phase", factory=MeshDistributor
    )
    distribute_mesh.meta["tip"] = "Distribute mesh among processes."

    insert_interfaces = pythia.pyre.inventory.facility(
        "insert_interfaces", family="initialize_phase", factory=MeshInsertInterfaces
    )
    insert_interfaces.meta["tip"] = (
        "Insert interfaces using PETSc mesh transform operation."
    )

    refine_mesh = pythia.pyre.inventory.facility(
        "refine_mesh", family="initialize_phase", factory=MeshRefiner
    )
    refine_mesh.meta["tip"] = (
        "Refine mesh."
    )

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="serial_phases"):
        """Constructor."""
        Component.__init__(self, name, facility="serial_phases")

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to guarantee order.

        """
        return [
            self.read_mesh,
            self.reorder_mesh,
            self.distribute_mesh,
            self.insert_interfaces,
            self.refine_mesh,
        ]



# End of file
