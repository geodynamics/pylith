#!/usr/bin/env nemesis
"""
Generate Gmsh files that contain common errors to ensure we trap them correctly.

1. Fault curve is not embedded in surface.
2. Fault surface is not embedded in domain.
3. Fault spline is not split at intersection.
"""

from abc import ABC, abstractmethod
import gmsh
import numpy

from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, group_exclude)

class GenerateApp(ABC):

    def __init__(self, cell, gui=False):
        self.cell = cell
        self.gui = gui
        self.geometry = None

    def main(self):
        self.create_geometry()
        self.mark()
        self.generate_mesh()
        self.write()
        self.finalize()

    @abstractmethod
    def create_geometry(self):
        pass

    @abstractmethod
    def mark(self):
        pass

    @abstractmethod
    def generate_mesh(self):
        pass

    def write(self):
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(f"{self.geometry}_{self.cell}.msh")

    def finalize(self):
        if self.gui:
            gmsh.fltk.run()
        gmsh.finalize()


class GenerateNoEmbed2D(GenerateApp):
    """
    p3-------p2
    |         |
    |  p4--p5 |
    |         |
    p0-------p1
    """

    DOMAIN_X = DOMAIN_Y = 8.0e+3
    DX = 2.0e+3

    def __init__(self, cell="tri", gui=False):
        super().__init__(cell, gui)
        self.geometry = "noembed"

    def create_geometry(self):

        gmsh.initialize()

        x0 = -0.5 * self.DOMAIN_X
        y0 = -0.5 * self.DOMAIN_Y
        p0 = gmsh.model.geo.add_point(x0, y0, 0)
        p1 = gmsh.model.geo.add_point(x0+self.DOMAIN_X, y0, 0)
        p2 = gmsh.model.geo.add_point(x0+self.DOMAIN_X, y0+self.DOMAIN_Y, 0)
        p3 = gmsh.model.geo.add_point(x0, y0+self.DOMAIN_Y, 0)
        p4 = gmsh.model.geo.add_point(x0+0.25*self.DOMAIN_X, y0+0.4*self.DOMAIN_Y, 0)
        p5 = gmsh.model.geo.add_point(x0+0.75*self.DOMAIN_X, y0+0.6*self.DOMAIN_Y, 0)

        self.c_yneg = gmsh.model.geo.add_line(p0, p1)
        self.c_xpos = gmsh.model.geo.add_line(p1, p2)
        self.c_ypos = gmsh.model.geo.add_line(p2, p3)
        self.c_xneg = gmsh.model.geo.add_line(p3, p0)
        self.c_fault = gmsh.model.geo.add_line(p4, p5)

        loop = gmsh.model.geo.add_curve_loop([self.c_yneg, self.c_xpos, self.c_ypos, self.c_xneg])
        self.surface = gmsh.model.geo.add_plane_surface([loop])

        gmsh.model.geo.synchronize()

    def mark(self):
        materials = (
            MaterialGroup(tag=1, entities=[self.surface]),
        )
        for material in materials:
            material.create_physical_group()

        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=10, dim=1, entities=[self.c_xneg]),
            VertexGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.c_xpos]),
            VertexGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.c_yneg]),
            VertexGroup(name="boundary_ypos", tag=13, dim=1, entities=[self.c_ypos]),
            VertexGroup(name="fault", tag=20, dim=1, entities=[self.c_fault]),
        )
        for group in vertex_groups:
            group.create_physical_group()

    def generate_mesh(self):
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        if self.cell == "quad":
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)

        gmsh.model.mesh.generate(2)


class GenerateNoEmbed3D(GenerateApp):

    DOMAIN_X = DOMAIN_Y = DOMAIN_Z = 8.0e+3
    DX = 2.0e+3

    def __init__(self, cell="tet", gui=False):
        super().__init__(cell, gui)
        self.geometry = "noembed"

    def create_geometry(self):

        gmsh.initialize()
        self.v_domain = gmsh.model.occ.add_box(-0.5*self.DOMAIN_X, -0.5*self.DOMAIN_Y, -self.DOMAIN_Z,
                                    self.DOMAIN_X, self.DOMAIN_Y, self.DOMAIN_Z)

        p0 = gmsh.model.occ.add_point(0.0, -0.2*self.DOMAIN_Y, 0)
        p1 = gmsh.model.occ.add_point(0.0, +0.2*self.DOMAIN_Y, 0)
        p2 = gmsh.model.occ.add_point(-0.1*self.DOMAIN_X, +0.2*self.DOMAIN_Y, -0.3*self.DOMAIN_Z)
        p3 = gmsh.model.occ.add_point(-0.1*self.DOMAIN_X, -0.2*self.DOMAIN_Y, -0.3*self.DOMAIN_Z)
        c_top = gmsh.model.occ.add_line(p0, p1)
        c_ypos = gmsh.model.occ.add_line(p1, p2)
        c_bot = gmsh.model.occ.add_line(p2, p3)
        c_yneg = gmsh.model.occ.add_line(p3, p0)
        loop = gmsh.model.occ.add_curve_loop([c_top, c_ypos, c_bot, c_yneg])
        self.s_fault = gmsh.model.occ.add_plane_surface([loop])

        gmsh.model.occ.synchronize()

    def mark(self):
        materials = (
            MaterialGroup(tag=1, entities=[self.v_domain]),
        )
        for material in materials:
            material.create_physical_group()

        vertex_groups = (
            VertexGroup(name="fault", tag=20, dim=2, entities=[self.s_fault]),
        )
        for group in vertex_groups:
            group.create_physical_group()

    def generate_mesh(self):
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.DX)
        if self.cell == "hex":
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)

        gmsh.model.mesh.generate(3)


class GenerateNoSplit2D(GenerateApp):
    """
    p3----------p2
    |            |
    |  p4-----p5 |
    |     p6     |
    |      |     |
    |     p7     |
    |            |
    p0----------p1
    """

    DOMAIN_X = DOMAIN_Y = 8.0e+3
    DX = 2.0e+3

    def __init__(self, cell="tri", gui=False):
        super().__init__(cell, gui)
        self.geometry = "nosplit"

    def create_geometry(self):

        gmsh.initialize()

        x0 = -0.5 * self.DOMAIN_X
        y0 = -0.5 * self.DOMAIN_Y
        p0 = gmsh.model.geo.add_point(x0, y0, 0)
        p1 = gmsh.model.geo.add_point(x0+self.DOMAIN_X, y0, 0)
        p2 = gmsh.model.geo.add_point(x0+self.DOMAIN_X, y0+self.DOMAIN_Y, 0)
        p3 = gmsh.model.geo.add_point(x0, y0+self.DOMAIN_Y, 0)

        self.c_yneg = gmsh.model.geo.add_line(p0, p1)
        self.c_xpos = gmsh.model.geo.add_line(p1, p2)
        self.c_ypos = gmsh.model.geo.add_line(p2, p3)
        self.c_xneg = gmsh.model.geo.add_line(p3, p0)
        loop = gmsh.model.geo.add_curve_loop([self.c_yneg, self.c_xpos, self.c_ypos, self.c_xneg])
        self.s_domain = gmsh.model.geo.add_plane_surface([loop])

        p4 = gmsh.model.geo.add_point(x0+0.25*self.DOMAIN_X, y0+0.4*self.DOMAIN_Y, 0)
        p5 = gmsh.model.geo.add_point(x0+0.75*self.DOMAIN_X, y0+0.6*self.DOMAIN_Y, 0)
        p6 = gmsh.model.geo.add_point(x0+0.4*self.DOMAIN_X, y0+0.52*self.DOMAIN_Y, 0)
        p7 = gmsh.model.geo.add_point(x0+0.54*self.DOMAIN_X, y0+0.51*self.DOMAIN_Y, 0)
        p8 = gmsh.model.geo.add_point(x0+0.4*self.DOMAIN_X, y0+0.3*self.DOMAIN_Y, 0)
        self.c_main_fault = gmsh.model.geo.add_spline([p4, p6, p7, p5])
        self.c_splay_fault = gmsh.model.geo.add_line(p7, p8)

        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(1, [self.c_main_fault], 2, self.s_domain)
        gmsh.model.mesh.embed(1, [self.c_splay_fault], 2, self.s_domain)

    def mark(self):
        materials = (
            MaterialGroup(tag=1, entities=[self.s_domain]),
        )
        for material in materials:
            material.create_physical_group()

        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=10, dim=1, entities=[self.c_xneg]),
            VertexGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.c_xpos]),
            VertexGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.c_yneg]),
            VertexGroup(name="boundary_ypos", tag=13, dim=1, entities=[self.c_ypos]),
            VertexGroup(name="main_fault", tag=20, dim=1, entities=[self.c_main_fault]),
            VertexGroup(name="splay_fault", tag=21, dim=1, entities=[self.c_splay_fault]),
        )
        for group in vertex_groups:
            group.create_physical_group()

    def generate_mesh(self):
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        if self.cell == "quad":
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)

        gmsh.model.mesh.generate(2)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", action="store", dest="geometry", required=True, choices=["noembed-2d","noembed-3d","nosplit-2d"])
    parser.add_argument("--cell", action="store", dest="cell", required=True, choices=["tri","quad","hex","tet"])
    parser.add_argument("--gui", action="store_true", dest="gui")
    args = parser.parse_args()

    if args.geometry == "noembed-2d":
        app = GenerateNoEmbed2D
    elif args.geometry == "noembed-3d":
        app = GenerateNoEmbed3D
    elif args.geometry == "nosplit-2d":
        app = GenerateNoSplit2D
    app(args.cell, args.gui).main()
