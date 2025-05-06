#!/usr/bin/env nemesis

import gmsh
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Block is DOMAIN_X by DOMAIN_Y with discretization size DX.

    p4----------p3
    |            |
    |            |
    |            |
    |            |
    p1----------p2
    """
    DOMAIN_X = DOMAIN_Y = 8.0e+3
    DX = 4.0e+3

    def __init__(self):
        self.cell_choices = {
            "required": True,
            "choices": ["tri", "quad"],
            }
        self.filename = "tri.msh"

    def create_geometry(self):
        """Create geometry.
        """
        lx = self.DOMAIN_X
        ly = self.DOMAIN_Y
        x0 = -0.5 * lx
        y0 = -0.5 * ly

        p1 = gmsh.model.geo.add_point(x0, y0, 0.0)
        p2 = gmsh.model.geo.add_point(x0+lx, y0, 0.0)
        p3 = gmsh.model.geo.add_point(x0+lx, y0+ly, 0.0)
        p4 = gmsh.model.geo.add_point(x0, y0+ly, 0.0)

        self.l_yneg = gmsh.model.geo.add_line(p1, p2)
        self.l_xpos = gmsh.model.geo.add_line(p2, p3)
        self.l_ypos = gmsh.model.geo.add_line(p3, p4)
        self.l_xneg = gmsh.model.geo.add_line(p4, p1)

        c0 = gmsh.model.geo.add_curve_loop([self.l_yneg, self.l_xpos, self.l_ypos, self.l_xneg])
        self.s_domain = gmsh.model.geo.add_plane_surface([c0])

        gmsh.model.geo.synchronize()

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.
        """
        materials = (
            MaterialGroup(tag=24, entities=[self.s_domain]),
        )
        for material in materials:
            material.create_physical_group()

        boundary_groups = (
            BoundaryGroup(name="boundary", tag=1, dim=1, entities=[self.l_xneg, self.l_yneg, self.l_xpos, self.l_ypos]),
        )
        for group in boundary_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh. Should also include optimizing the mesh quality.
        """
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.DX)
        if cell == "quad":
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)

        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()


# End of file
