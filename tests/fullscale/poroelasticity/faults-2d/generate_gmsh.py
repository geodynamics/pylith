#!/usr/bin/env nemesis

import gmsh
from pylith.meshio.gmsh_utils import BoundaryGroup, MaterialGroup, GenerateMesh


class App(GenerateMesh):
    """
    Block is DOMAIN_X by DOMAIN_Y with discretization size DX.

    p4------------p3
    |              |
    |              |
    p5------------p6
    |              |
    |              |
    p1------------p2
    """

    DOMAIN_X = DOMAIN_Y = 8.0e3
    DX = 2.0e3

    def __init__(self):
        self.cell_choices = {
            "required": True,
            "choices": ["tri", "quad"],
        }
        self.filename = "mesh.msh"

    def create_geometry(self):
        """Create geometry."""
        lx = self.DOMAIN_X
        ly = self.DOMAIN_Y
        x0 = 0.0
        y0 = -0.5 * ly

        p1 = gmsh.model.geo.add_point(x0, y0, 0.0)
        p2 = gmsh.model.geo.add_point(x0 + lx, y0, 0.0)
        p3 = gmsh.model.geo.add_point(x0 + lx, y0 + ly, 0.0)
        p4 = gmsh.model.geo.add_point(x0, y0 + ly, 0.0)

        p5 = gmsh.model.geo.add_point(x0, y0 + 0.5 * ly, 0.0)
        p6 = gmsh.model.geo.add_point(x0 + lx, y0 + 0.5 * ly, 0.0)

        self.xmid_corners = [p5, p6]

        self.l_yneg = gmsh.model.geo.add_line(p1, p2)
        self.l_xpos0 = gmsh.model.geo.add_line(p2, p6)
        self.l_xpos1 = gmsh.model.geo.add_line(p6, p3)
        self.l_ypos = gmsh.model.geo.add_line(p3, p4)
        self.l_xneg0 = gmsh.model.geo.add_line(p4, p5)
        self.l_xneg1 = gmsh.model.geo.add_line(p5, p1)
        self.l_fault = gmsh.model.geo.add_line(p5, p6)

        c1 = gmsh.model.geo.add_curve_loop(
            [self.l_yneg, self.l_xpos0, -self.l_fault, self.l_xneg1]
        )
        self.s_xneg = gmsh.model.geo.add_plane_surface([c1])
        c2 = gmsh.model.geo.add_curve_loop(
            [self.l_fault, self.l_xpos1, self.l_ypos, self.l_xneg0]
        )
        self.s_xpos = gmsh.model.geo.add_plane_surface([c2])

        gmsh.model.geo.synchronize()

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc."""
        materials = (MaterialGroup(tag=1, entities=[self.s_xneg, self.s_xpos]),)
        for material in materials:
            material.create_physical_group()

        face_groups = (
            BoundaryGroup(
                name="boundary_xneg",
                tag=10,
                dim=1,
                entities=[self.l_xneg0, self.l_xneg1],
            ),
            BoundaryGroup(
                name="boundary_xpos",
                tag=11,
                dim=1,
                entities=[self.l_xpos0, self.l_xpos1],
            ),
            BoundaryGroup(
                name="boundary_yneg",
                tag=12,
                dim=1,
                entities=[self.l_yneg],
            ),
            BoundaryGroup(
                name="boundary_ypos",
                tag=13,
                dim=1,
                entities=[self.l_ypos],
            ),
            BoundaryGroup(name="fault", tag=20, dim=1, entities=[self.l_fault]),
        )
        for group in face_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh. Should also include optimizing the mesh quality."""
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.DX)
        if cell == "quad":
            # Generate a tri mesh and then recombine cells to form quadrilaterals.
            # We use the Frontal-Delaunay for Quads algorithm.
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.recombine()
        else:
            gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()


# End of file
