#!/usr/bin/env nemesis

import numpy
import gmsh
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Block is DOMAIN_X by DOMAIN_Y x DOMAN_Z with discretization size DX.

    -30km <= x <= +30km
    -30km <= y <= +30km
    -30km <= z <= 0
    """
    DOMAIN_X = DOMAIN_Y = 40.0e+3
    DOMAIN_Z = 20.0e+3
    FAULT_LENGTH = 20.0e+3
    FAULT_WIDTH = 10.0e+3
    DX = 5.0e+3

    def __init__(self):
        self.cell_choices = {
            "required": False,
            "default": "tet",
            "choices": ["tet"],
            }
        self.filename = "mesh_tet.msh"

    def create_geometry(self):
        """Create geometry.
        """
        gmsh.initialize()
        self.box = gmsh.model.occ.add_box(-0.5*self.DOMAIN_X, -0.5*self.DOMAIN_Y, -self.DOMAIN_Z,
                                    self.DOMAIN_X, self.DOMAIN_Y, self.DOMAIN_Z)

        self.s_xneg = 1
        self.s_xpos = 2
        self.s_yneg = 3
        self.s_ypos = 4
        self.s_zneg = 5
        self.s_zpos = 6
        
        p_faulttrace_yneg = gmsh.model.occ.add_point(0.0, -0.5*self.FAULT_LENGTH, 0.0)
        p_faulttrace_ypos = gmsh.model.occ.add_point(0.0, +0.5*self.FAULT_LENGTH, 0.0)
        p_faultbot_yneg = gmsh.model.occ.add_point(0.0, -0.5*self.FAULT_LENGTH, -self.FAULT_WIDTH)
        p_faultbot_ypos = gmsh.model.occ.add_point(0.0, +0.5*self.FAULT_LENGTH, -self.FAULT_WIDTH)

        self.c_faulttrace = gmsh.model.occ.add_line(p_faulttrace_yneg, p_faulttrace_ypos)
        self.c_fault_yneg = gmsh.model.occ.add_line(p_faulttrace_yneg, p_faultbot_yneg)
        self.c_faultbot = gmsh.model.occ.add_line(p_faultbot_yneg, p_faultbot_ypos)
        self.c_fault_ypos = gmsh.model.occ.add_line(p_faulttrace_ypos, p_faultbot_ypos)

        loop_fault = gmsh.model.occ.add_curve_loop([-self.c_faulttrace, self.c_fault_yneg, self.c_faultbot, -self.c_fault_ypos])
        self.s_fault = gmsh.model.occ.add_plane_surface([loop_fault])
        
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.embed(2, [self.s_fault], 3, self.box)
        gmsh.model.mesh.embed(1, [self.c_faulttrace], 2, self.s_zpos)

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.
        """
        materials = (
            MaterialGroup(tag=1, entities=[self.box]),
        )
        for material in materials:
            material.create_physical_group()

        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=10, dim=2, entities=[self.s_xneg]),
            VertexGroup(name="boundary_xpos", tag=11, dim=2, entities=[self.s_xpos]),
            VertexGroup(name="boundary_yneg", tag=12, dim=2, entities=[self.s_yneg]),
            VertexGroup(name="boundary_ypos", tag=13, dim=2, entities=[self.s_ypos]),
            VertexGroup(name="boundary_zneg", tag=14, dim=2, entities=[self.s_zneg]),
            VertexGroup(name="boundary_zpos", tag=15, dim=2, entities=[self.s_zpos]),
            VertexGroup(name="fault", tag=20, dim=2, entities=[self.s_fault]),
            VertexGroup(name="fault_edges", tag=21, dim=1, entities=[self.c_fault_yneg, self.c_faultbot, self.c_fault_ypos]),
        )
        for group in vertex_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh. Should also include optimizing the mesh quality.
        """
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.DX)
        if cell == "hex":
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
            #nnodes = 3 + int(self.DOMAIN_Y/self.DX + 0.5)
            #gmsh.model.mesh.set_transfinite_curve(self.l_xneg, nnodes)
            #gmsh.model.mesh.set_transfinite_curve(self.l_split0, nnodes)
            #gmsh.model.mesh.set_transfinite_surface(self.s_xmid, cornerTags=self.xmid_corners)
            #gmsh.model.mesh.set_recombine(3, self.s_xmid)

        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()


# End of file
