#!/usr/bin/env python3

import numpy
import gmsh
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Block is DOMAIN_X by DOMAIN_Y x DOMAN_Z with discretization size DX.

    -4km <= x <= +4km
    -4km <= y <= +4km
    -8km <= z <= 0
    """
    DOMAIN_X = DOMAIN_Y = DOMAIN_Z = 32.0e+3
    DX = 4.0e+3

    def __init__(self):
        self.cell_choices = {
            "required": True,
            "choices": ["tet", "hex"],
            }
        self.filename = "mesh.msh"

    def create_geometry(self):
        """Create geometry.
        """
        gmsh.initialize()
        gmsh.model.add("mesh")
        box = gmsh.model.occ.add_box(-0.5*self.DOMAIN_X, -0.5*self.DOMAIN_Y, -self.DOMAIN_Z,
                                    self.DOMAIN_X, self.DOMAIN_Y, self.DOMAIN_Z)

        disk = gmsh.model.occ.add_disk(0.0, 0.0, -0.5*self.DOMAIN_Z, self.DOMAIN_Y, self.DOMAIN_Z)
        gmsh.model.occ.rotate([(2,disk)], 0.0, 0.0, -0.5*self.DOMAIN_Z, 0.0, 1.0, 0.0, 0.4*numpy.pi)
        dimTags, dimTagsMap = gmsh.model.occ.fragment([(3,box)], [(2,disk)])
        gmsh.model.occ.remove([dimTags[-1]], recursive=True)

        gmsh.model.occ.synchronize()

        bbox = gmsh.model.get_bounding_box(-1, -1)
        dx = 100.0
        self.s_xneg = gmsh.model.get_entities_in_bounding_box(bbox[0]-dx, bbox[1]-dx, bbox[2]-dx, bbox[0]+dx, bbox[4]+dx, bbox[5]+dx, dim=2)
        self.s_xpos = gmsh.model.get_entities_in_bounding_box(bbox[3]-dx, bbox[1]-dx, bbox[2]-dx, bbox[3]+dx, bbox[4]+dx, bbox[5]+dx, dim=2)
        self.s_yneg = gmsh.model.get_entities_in_bounding_box(bbox[0]-dx, bbox[1]-dx, bbox[2]-dx, bbox[3]+dx, bbox[1]+dx, bbox[5]+dx, dim=2)
        self.s_ypos = gmsh.model.get_entities_in_bounding_box(bbox[0]-dx, bbox[4]-dx, bbox[2]-dx, bbox[3]+dx, bbox[4]+dx, bbox[5]+dx, dim=2)
        self.s_zneg = gmsh.model.get_entities_in_bounding_box(bbox[0]-dx, bbox[1]-dx, bbox[2]-dx, bbox[3]+dx, bbox[4]+dx, bbox[2]+dx, dim=2)
        self.s_zpos = gmsh.model.get_entities_in_bounding_box(bbox[0]-dx, bbox[1]-dx, bbox[5]-dx, bbox[3]+dx, bbox[4]+dx, bbox[5]+dx, dim=2)
        xtol = 0.2*self.DOMAIN_X
        self.s_fault = gmsh.model.get_entities_in_bounding_box(-xtol, bbox[1]-dx, bbox[2]-dx, +xtol, bbox[4]+dx, bbox[5]+dx, dim=2)

        self.v_xneg = dimTags[0]
        self.v_xpos = dimTags[1]

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.
        """
        materials = (
            MaterialGroup(tag=1, entities=[tag for dim, tag in [self.v_xneg, self.v_xpos]]),
        )
        for material in materials:
            material.create_physical_group()

        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=10, dim=2, entities=[tag for dim, tag in self.s_xneg]),
            VertexGroup(name="boundary_xpos", tag=11, dim=2, entities=[tag for dim, tag in self.s_xpos]),
            VertexGroup(name="boundary_yneg", tag=12, dim=2, entities=[tag for dim, tag in self.s_yneg]),
            VertexGroup(name="boundary_ypos", tag=13, dim=2, entities=[tag for dim, tag in self.s_ypos]),
            VertexGroup(name="boundary_zneg", tag=14, dim=2, entities=[tag for dim, tag in self.s_zneg]),
            VertexGroup(name="boundary_zpos", tag=15, dim=2, entities=[tag for dim, tag in self.s_zpos]),
            VertexGroup(name="fault", tag=20, dim=2, entities=[tag for dim,tag in self.s_fault]),
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
