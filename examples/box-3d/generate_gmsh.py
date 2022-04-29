#!/usr/bin/env nemesis
"""Generate a hex or tet mesh using Gmsh, making use of the OpenCascade geometry engine.

We use the `gmsh_utils` module provided with PyLith. This module has helper functions
and classes for marking materials and boundaries compatible with PyLith. We also use the
`GenerateMesh` class from the module as a base class to our local `App` class. The
`GenerateMesh` class handles processing of command line options, initialization, and
finalizing the mesh.

Run `generate_gmsh.py --help` to see the command line options.
"""

import gmsh
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Application for generating the mesh.

    Block is 12km x 12km x 9km
    -6 km <= x <= 6 km
    -6 km <= y <= 6 km
    -9 km <= z <= 0 km

    Discretization size is 3km.
    """
    DOMAIN_X = DOMAIN_Y = 12.0e+3
    DOMAIN_Z = 9.0e+3
    DX = 3.0e+3

    def __init__(self):
        """Constructor.
        """
        # Set the cell choices available through command line options
        # with the default cell type `hex` matching the PyLith parameter files.
        self.cell_choices = {
            "default": "hex",
            "choices": ["tet", "hex"],
            }

    def create_geometry(self):
        """Create geometry.

        This method is abstract in the base class and must be implemented.
        """
        # Define the coordinates of one corner.
        x0 = -0.5 * self.DOMAIN_X
        y0 = -0.5 * self.DOMAIN_Y
        z0 = -self.DOMAIN_Z

        # Create a box
        self.v1 = gmsh.model.occ.add_box(x0, y0, z0, self.DOMAIN_X, self.DOMAIN_Y, self.DOMAIN_Z)

        gmsh.model.occ.synchronize()

        # Get bounding surfaces using bounding boxes. This makes the script completely independent
        # of how the surfaces are numbered, which can change.
        dx = 100.0
        bbox = gmsh.model.get_bounding_box(-1, -1)
        get_bbox_entities = gmsh.model.get_entities_in_bounding_box
        self.s_xneg = get_bbox_entities(bbox[0]-dx, bbox[1]-dx, bbox[2]-dx, bbox[0]+dx, bbox[4]+dx, bbox[5]+dx, dim=2)
        self.s_xpos = get_bbox_entities(bbox[3]-dx, bbox[1]-dx, bbox[2]-dx, bbox[3]+dx, bbox[4]+dx, bbox[5]+dx, dim=2)
        self.s_yneg = get_bbox_entities(bbox[0]-dx, bbox[1]-dx, bbox[2]-dx, bbox[3]+dx, bbox[1]+dx, bbox[5]+dx, dim=2)
        self.s_ypos = get_bbox_entities(bbox[0]-dx, bbox[4]-dx, bbox[2]-dx, bbox[3]+dx, bbox[4]+dx, bbox[5]+dx, dim=2)
        self.s_zneg = get_bbox_entities(bbox[0]-dx, bbox[1]-dx, bbox[2]-dx, bbox[3]+dx, bbox[4]+dx, bbox[2]+dx, dim=2)
        self.s_zpos = get_bbox_entities(bbox[0]-dx, bbox[1]-dx, bbox[5]-dx, bbox[3]+dx, bbox[4]+dx, bbox[5]+dx, dim=2)

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented.
        """
        # We use the `MaterialGroup` data class and `create_material` function from `gmsh_utils`
        # to define the materials. The `MaterialGroup` constructor takes the tag of the physical
        # group to create for the material as its first argument and a list of geometry entities
        # (surface in 2D and volume in 3D) as its second argument. The tag of the physical 
        # group corresponds to the `label_value` Pyre property for the material in the PyLith
        # settings.
        materials = (
            MaterialGroup(tag=1, entities=[self.v1]),
        )
        for material in materials:
            material.create_physical_group()

        # We use the `VertexGroup` data class and `create_group` function from `gmsh_utils`
        # to define the groups of vertices. The `VertexGroup` constructor takes the name and
        # tag of the physical group as the first two arguments and the list of entities
        # for the group as the third argument. The name corresponds to the `label` Pyre property
        # and the tag corresponds to the `label_value` Pyre property in the boundary condition
        # settings.
        #
        # The `create_group` function will automatically create the physical groups for all
        # lower dimension entities; these are needed by PyLith.
        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=10, dim=2, entities=[tag for dim, tag in self.s_xneg]),
            VertexGroup(name="boundary_xpos", tag=11, dim=2, entities=[tag for dim, tag in self.s_xpos]),
            VertexGroup(name="boundary_yneg", tag=12, dim=2, entities=[tag for dim, tag in self.s_yneg]),
            VertexGroup(name="boundary_ypos", tag=13, dim=2, entities=[tag for dim, tag in self.s_ypos]),
            VertexGroup(name="boundary_zneg", tag=14, dim=2, entities=[tag for dim, tag in self.s_zneg]),
            VertexGroup(name="boundary_zpos", tag=15, dim=2, entities=[tag for dim, tag in self.s_zpos]),
        )
        for group in vertex_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh.

        This method is abstract in the base class and must be implemented.
        """
        # Set the cell size
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.DX)
        if cell == "hex":
            # A transfinite mesh with recombine set to True should result in a hex mesh. 
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)

        # Generate the mesh and then improve mesh quality using Laplacian smoothing.
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()


# End of file
