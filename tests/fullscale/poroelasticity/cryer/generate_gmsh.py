#!/usr/bin/env nemesis

import math

import gmsh
from pylith.meshio.gmsh_utils import BoundaryGroup, MaterialGroup, GenerateMesh


class App(GenerateMesh):
    """
    Domain is a sphere with radius RADIUS and finest discretization size DX.
    """

    RADIUS = 10.0e+3
    DX = 1.0e+3
    
    def __init__(self):
        self.cell_choices = {
            "required": False,
            "choices": ["tet"],
        }
        self.filename = "mesh_tet.msh"

    def create_geometry(self):
        """Create geometry."""
        x0 = 0.0
        y0 = 0.0
        z0 = 0.0

        self.v_domain = gmsh.model.occ.add_sphere(x0, y0, z0, radius=self.RADIUS, angle1=0, angle3=math.pi/2)
        self.s_shell = 1
        self.s_zneg = 2
        self.s_yneg = 3
        self.s_xneg = 4
        gmsh.model.occ.synchronize()

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc."""
        materials = (MaterialGroup(tag=1, entities=[self.v_domain]),)
        for material in materials:
            material.create_physical_group()

        face_groups = (
            BoundaryGroup(
                name="boundary_shell",
                tag=10,
                dim=2,
                entities=[self.s_shell],
            ),
            BoundaryGroup(
                name="boundary_shell_copy",
                tag=11,
                dim=2,
                entities=[self.s_shell],
            ),
            BoundaryGroup(
                name="boundary_zneg",
                tag=12,
                dim=2,
                entities=[self.s_zneg],
            ),
            BoundaryGroup(
                name="boundary_yneg",
                tag=13,
                dim=2,
                entities=[self.s_yneg],
            ),
            BoundaryGroup(
                name="boundary_xneg",
                tag=14,
                dim=2,
                entities=[self.s_xneg],
            ),
        )
        for group in face_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh. Should also include optimizing the mesh quality."""
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.DX)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()
