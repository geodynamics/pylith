#!/usr/bin/env nemesis
"""
Physical groups have a name, dimension, and physical tag. We will make a label for each physical group.
Entities have a dimension, bounding box, and physical tag
Nodes have a node tag, coordinates, and an entity. We need to lookup the entity to get the physical tag for a node.
Elements have an element tag, dimension, element type, list of node tags, and an entity. We need to lookup the entity to get the physical tag for a node.

v5----v4----v3
|      |     |
|      |     |
|      |     |
v0----v1----v2

"""

import gmsh
import numpy

def create_group(name, dim, entities, recursive=True, tag=-1):
    tag = gmsh.model.addPhysicalGroup(dim, entities, tag)
    gmsh.model.setPhysicalName(dim, tag, name)
    entities_lowerdim = []
    for entity in entities:
        entities_up, entities_down = gmsh.model.get_adjacencies(dim, entity)
        entities_lowerdim += [e for e in entities_down]
    if recursive and dim >= 1:
        create_group(name, dim-1, entities_lowerdim, tag=tag)

def create_material(tag, entities):
    if tag <= 0:
        raise ValueError(f"ERROR: Attempting to use non-positive material tag '{tag}'. Tags for physical groups must be positive.")
    dim = gmsh.model.get_dimension()
    name = gmsh.model.get_physical_name(dim, tag)
    if name:
        raise ValueError(f"ERROR: Attempting to use material tag '{tag}' that is already in use for material '{name}'.")
    gmsh.model.addPhysicalGroup(dim, entities, tag)
    gmsh.model.setPhysicalName(dim, tag, f"material-id:{tag}")

class Generate2D():

    DOMAIN_X = DOMAIN_Y = 8.0e+3
    DX = 1.0e+3

    def __init__(self, cell="tri"):
        self.cell = cell

    def main(self):
        self.create_geometry()
        self.mark()
        self.generate_mesh()
        self.write()
        self.finalize()

    def create_geometry(self):

        gmsh.initialize()
        gmsh.model.add("mesh")

        x0 = -0.5 * self.DOMAIN_X
        y0 = -0.5 * self.DOMAIN_Y
        p0 = gmsh.model.geo.addPoint(x0, y0, 0, self.DX)
        p1 = gmsh.model.geo.addPoint(x0+0.5*self.DOMAIN_X, y0, 0, self.DX)
        p2 = gmsh.model.geo.addPoint(x0+self.DOMAIN_X, y0, 0, self.DX)
        p3 = gmsh.model.geo.addPoint(x0+self.DOMAIN_X, y0+self.DOMAIN_Y, 0, self.DX)
        p4 = gmsh.model.geo.addPoint(x0+0.5*self.DOMAIN_X, y0+self.DOMAIN_Y, 0, self.DX)
        p5 = gmsh.model.geo.addPoint(x0, y0+self.DOMAIN_Y, 0, self.DX)

        self.l_yneg0 = gmsh.model.geo.addLine(p0, p1)
        self.l_yneg1 = gmsh.model.geo.addLine(p1, p2)
        self.l_xpos = gmsh.model.geo.addLine(p2, p3)
        self.l_ypos0 = gmsh.model.geo.addLine(p3, p4)
        self.l_ypos1 = gmsh.model.geo.addLine(p4, p5)
        self.l_xneg = gmsh.model.geo.addLine(p5, p0)
        self.l_fault = gmsh.model.geo.addLine(p1, p4)

        c0 = gmsh.model.geo.addCurveLoop([self.l_yneg0, self.l_fault, self.l_ypos1, self.l_xneg])
        self.s_xneg = gmsh.model.geo.addPlaneSurface([c0])
        c1 = gmsh.model.geo.addCurveLoop([self.l_yneg1, self.l_xpos, self.l_ypos0, -self.l_fault])
        self.s_xpos = gmsh.model.geo.addPlaneSurface([c1])

        gmsh.model.geo.synchronize()

    def mark(self):
        materials = (
            (1, (self.s_xneg,)),
            (2, (self.s_xpos,)),
        )
        for matid, material in materials:
            create_material(matid, material)

        edge_groups = {
            "fault": (self.l_fault,),
            "boundary_xneg": (self.l_xneg,),
            "boundary_xpos": (self.l_xpos,),
            "boundary_yneg": (self.l_yneg0, self.l_yneg1),
            "boundary_ypos": (self.l_ypos0, self.l_ypos1),
        }
        for name, group in edge_groups.items():
            create_group(name, gmsh.model.get_dimension()-1, group)

    def generate_mesh(self):
        # right angle tris
        gmsh.option.setNumber("Mesh.Algorithm", 8)

        if self.cell == "quad":
            gmsh.option.setNumber("Mesh.RecombineAll", 1)

        # Can set refinement criterion
        # gmsh.model.mesh.setSize(0, maxSize)

        gmsh.model.mesh.generate(2)

    def write(self):
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(f"{self.cell}_ascii.msh")

        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(f"{self.cell}_binary.msh")

    def finalize(self):
        gmsh.fltk.run()
        gmsh.finalize()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--cell", action="store", dest="cell", required=True, choices=["tri","quad"])
    args = parser.parse_args()

    Generate2D(args.cell).main()
