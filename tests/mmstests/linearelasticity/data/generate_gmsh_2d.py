#!/usr/bin/env python3

import gmsh

DOMAIN_X = DOMAIN_Y = 8.0e+3
DX = 4.0e+3

# Physical groups have a name, dimension, and physical tag. We will make a label for each physical group.
# Entities have a dimension, bounding box, and physical tag
# Nodes have a node tag, coordinates, and an entity. We need to lookup the entity to get the physical tag for a node.
# Elements have an element tag, dimension, element type, list of node tags, and an entity. We need to lookup the entity to get the physical tag for a node.

def generate_mesh(cell="tri"):

    gmsh.initialize()
    gmsh.model.add("mesh")

    x0 = -0.5 * DOMAIN_X
    y0 = -0.5 * DOMAIN_Y
    v0 = gmsh.model.geo.addPoint(x0, y0, 0, DX)
    v1 = gmsh.model.geo.addPoint(x0+DOMAIN_X, y0, 0, DX)
    v2 = gmsh.model.geo.addPoint(x0+DOMAIN_X, y0+DOMAIN_Y, 0, DX)
    v3 = gmsh.model.geo.addPoint(x0, y0+DOMAIN_Y, 0, DX)

    e0 = gmsh.model.geo.addLine(v0, v1)
    e1 = gmsh.model.geo.addLine(v1, v2)
    e2 = gmsh.model.geo.addLine(v2, v3)
    e3 = gmsh.model.geo.addLine(v3, v0)

    c0 = gmsh.model.geo.addCurveLoop([e0, e1, e2, e3])
    s0 = gmsh.model.geo.addPlaneSurface([c0])

    gmsh.model.geo.synchronize()

    # This adds "boundary" to the PhysicalNames section with the next available tag
    boundaryV = gmsh.model.addPhysicalGroup(0, [v0, v1, v2, v3], 1)
    gmsh.model.setPhysicalName(0, boundaryV, "boundary")

    # This adds "boundary" to the PhysicalNames section with the next available tag
    boundary = gmsh.model.addPhysicalGroup(1, [e0, e1, e2, e3], 1)
    gmsh.model.setPhysicalName(1, boundary, "boundary")

    # This adds "material-id" to the PhysicalNames section with tag 24
    material = gmsh.model.addPhysicalGroup(2, [s0], 24)
    gmsh.model.setPhysicalName(2, material, "material-id")

    # right angle tris
    gmsh.option.setNumber("Mesh.Algorithm", 8)

    if cell == "quad":
        gmsh.option.setNumber("Mesh.RecombineAll", 1)

    # Can set refinement criterion
    # gmsh.model.mesh.setSize(0, maxSize)

    gmsh.model.mesh.generate(2)
    gmsh.write(f"{cell}_gmsh.msh")

    gmsh.fltk.run()
    gmsh.finalize()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--cell", action="store", dest="cell", required=True, choices=["tri","quad"])
    args = parser.parse_args()

    generate_mesh(args.cell)
