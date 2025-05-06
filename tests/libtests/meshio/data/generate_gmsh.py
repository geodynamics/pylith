#!/usr/bin/env nemesis
"""
Physical groups have a name, dimension, and physical tag. We will make a label for each physical group.
Entities have a dimension, bounding box, and physical tag
Nodes have a node tag, coordinates, and an entity. We need to lookup the entity to get the physical tag for a node.
Elements have an element tag, dimension, element type, list of node tags, and an entity. We need to lookup the entity to get the physical tag for a node.
"""

from abc import ABC, abstractmethod
import gmsh
import numpy

from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup)


def get_vertices():
    dim = gmsh.model.get_dimension()
    node_tags, node_coords, node_params = gmsh.model.mesh.get_nodes()
    num_nodes = node_coords.shape[0] // 3
    node_coords = node_coords.reshape((num_nodes, 3))
    return node_coords[:,0:dim]


def get_cells():
    elem_types, elem_tags, node_tags = gmsh.model.mesh.get_elements()
    if elem_types[-2] == 2:
        ncorners = 3
    elif elem_types[-2] == 3:
        ncorners = 4
    elif elem_types[-2] == 4:
        ncorners = 4
    elif elem_types[-2] == 5:
        ncorners = 8
    else:
        raise ValueError(f"Unknown cell type '{elem_types[-2]}'.")
    num_cells = node_tags[-2].shape[0] // ncorners
    cells = node_tags[-2].reshape((num_cells, ncorners))
    return cells


def get_material_ids(cells):
    dim = gmsh.model.get_dimension()
    physical_groups = gmsh.model.get_physical_groups(dim)
    cell_indices = []
    min_tag = numpy.uint64(-1)
    for dim, tag in physical_groups:
        group_entities = gmsh.model.get_entities_for_physical_group(dim, tag)
        indices_tag = []
        for group_entity in group_entities:
            entities = gmsh.model.mesh.get_elements(dim, group_entity)
            indices_tag += entities[1][0].tolist()
            min_tag = min(numpy.min(entities[1][0]), min_tag)
        cell_indices.append(indices_tag)
    material_ids = numpy.zeros(cells.shape[0], dtype=numpy.uint64)
    for (dim, tag), indices in zip(physical_groups, cell_indices):
        indices -= min_tag
        material_ids[indices] = tag
    return material_ids


def get_vertex_groups():
    """IMPORTANT: This grabs the nodes at the highest dimension. It does not produce
    correct results if we have excluded entities at lower dimensions."""
    mesh_dim = gmsh.model.get_dimension()
    groups = {}
    for dim in range(mesh_dim-1,-1,-1):
        physical_groups = gmsh.model.get_physical_groups(dim)
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            node_tags, node_coords = gmsh.model.mesh.get_nodes_for_physical_group(dim, tag)
            if not name in groups:
                groups[name] = node_tags
    return groups


class WriterCxx():

    def __init__(self, filename="stdout"):
        if filename == "stdout":
            import sys
            self.fout = sys.stdout
        else:
            self.fout = open(filename, "w")
    
    def write_vertices(self, vertices):
        nvertices, space_dim = vertices.shape
        self.fout.write(f"const size_t numVertices = {nvertices};\n")
        self.fout.write(f"const size_t space_dim = {space_dim};\n")
        self.fout.write("static const PylithScalar vertices[numVertices*spaceDim] = " + "{\n")
        for vertex in vertices:
            line = ", ".join([f"{x:+16.8e}" for x in vertex]) + ",\n"
            self.fout.write(line)
        self.fout.write("};\n")
        self.fout.write("delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);\n\n")

    def write_cells(self, cells):
        ncells, ncorners = cells.shape
        self.fout.write(f"const size_t numCells = {ncells};\n")
        self.fout.write(f"const size_t numCorners = {ncorners};\n")
        self.fout.write("static const PylithInt cells[numCells*numCorners] = " + "{\n")
        for cell in cells[:,0:ncorners]:
            line = ", ".join([f"{int(node-1):d}" for node in cell]) + ",\n"
            self.fout.write(line)
        self.fout.write("};\n")
        self.fout.write("delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);\n\n")

    def write_material_ids(self, material_ids):
        ncells = len(material_ids)
        self.fout.write("static const PylithInt materialIds[numCells] = " + "{\n")
        line = [f"{mid}," for mid in material_ids]
        self.fout.write("\n".join(line))
        self.fout.write("\n};\n")
        self.fout.write("_data->materialIds = const_cast<PylithInt*>(materialIds);\n\n")

    def write_vertex_groups(self, groups):
        ngroups = len(groups)
        self.fout.write(f"_data->numVertexGroups = {ngroups};\n")
        gsizes = [len(group) for group in groups.values()]
        gsizes_str = [f"{len(group)}" for group in groups.values()]
        self.fout.write(f"static const PylithInt vertexGroupSizes[{ngroups}] = ")
        self.fout.write("{" + ", ".join(gsizes_str) + "};\n")
        self.fout.write("_data->vertexGroupSizes = const_cast<PylithInt*>(vertexGroupSizes);\n")
        self.fout.write(f"static const PylithInt vertexGroups[{sum(gsizes)}] = " + "{\n")
        for group in groups.values():
            line = ", ".join([f"{int(vertex-1):d}" for vertex in group]) + ",\n"
            self.fout.write(line)
        self.fout.write("};\n")
        self.fout.write("_data->vertexGroups = const_cast<PylithInt*>(vertexGroups);\n")
        self.fout.write(f"static const char* vertexGroupNames[{ngroups}] = " + "{\n")
        for name in groups.keys():
            self.fout.write(f'"{name}",\n')
        self.fout.write("};\n")
        self.fout.write("_data->vertexGroupNames = const_cast<char**>(vertexGroupNames);\n")


def write_mesh(writer=WriterCxx()):
    writer.write_vertices(get_vertices())
    cells = get_cells()
    writer.write_cells(cells)
    writer.write_material_ids(get_material_ids(cells))
    del cells
    writer.write_vertex_groups(get_vertex_groups())


class GenerateApp(ABC):

    def __init__(self, cell, mark_vertices=False, dump_mesh=False, gui=False):
        self.cell = cell
        self.dump_mesh = dump_mesh
        self.gui = gui
        self.mark_vertices = mark_vertices
        self.geometry = None

    def main(self):
        self.create_geometry()
        self.mark(self.mark_vertices)
        self.generate_mesh()
        self.write()
        if self.dump_mesh:
            self.print_mesh()
        self.finalize()

    @abstractmethod
    def create_geometry(self):
        pass

    @abstractmethod
    def mark(self, mark_vertices):
        pass

    @abstractmethod
    def generate_mesh(self):
        pass

    def write(self):
        mark_flag = "vertices" if self.mark_vertices else "boundary"
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(f"{self.geometry}_{self.cell}_{mark_flag}_ascii.msh")

        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(f"{self.geometry}_{self.cell}_{mark_flag}_binary.msh")

    def print_mesh(self):
        write_mesh(writer=WriterCxx())

    def finalize(self):
        if self.gui:
            gmsh.fltk.run()
        gmsh.finalize()


class GenerateBox2D(GenerateApp):
    """
    v5----v4----v3
    |      |     |
    |      |     |
    |      |     |
    v0----v1----v2
    """

    DOMAIN_X = DOMAIN_Y = 8.0e+3
    DX = 4.0e+3

    def __init__(self, cell="tri", mark_vertices=False, dump_mesh=False, gui=False):
        super().__init__(cell, mark_vertices, dump_mesh, gui)
        self.geometry = "box"

    def create_geometry(self):

        gmsh.initialize()
        gmsh.model.add("mesh")

        x0 = -0.5 * self.DOMAIN_X
        y0 = -0.5 * self.DOMAIN_Y
        p0 = gmsh.model.geo.addPoint(x0, y0, 0)
        p1 = gmsh.model.geo.addPoint(x0+0.5*self.DOMAIN_X, y0, 0)
        p2 = gmsh.model.geo.addPoint(x0+self.DOMAIN_X, y0, 0)
        p3 = gmsh.model.geo.addPoint(x0+self.DOMAIN_X, y0+self.DOMAIN_Y, 0)
        p4 = gmsh.model.geo.addPoint(x0+0.5*self.DOMAIN_X, y0+self.DOMAIN_Y, 0)
        p5 = gmsh.model.geo.addPoint(x0, y0+self.DOMAIN_Y, 0)

        self.l_yneg0 = gmsh.model.geo.addLine(p0, p1)
        self.l_yneg1 = gmsh.model.geo.addLine(p1, p2)
        self.l_xpos = gmsh.model.geo.addLine(p2, p3)
        self.l_ypos0 = gmsh.model.geo.addLine(p3, p4)
        self.l_ypos1 = gmsh.model.geo.addLine(p4, p5)
        self.l_xneg = gmsh.model.geo.addLine(p5, p0)
        self.l_fault = gmsh.model.geo.addLine(p1, p4)
        self.p_fault_end = p1

        c0 = gmsh.model.geo.addCurveLoop([self.l_yneg0, self.l_fault, self.l_ypos1, self.l_xneg])
        self.s_xneg = gmsh.model.geo.addPlaneSurface([c0])
        c1 = gmsh.model.geo.addCurveLoop([self.l_yneg1, self.l_xpos, self.l_ypos0, -self.l_fault])
        self.s_xpos = gmsh.model.geo.addPlaneSurface([c1])

        gmsh.model.geo.synchronize()

    def mark(self, mark_vertices):
        materials = (
            MaterialGroup(tag=1, entities=[self.s_xpos]),
            MaterialGroup(tag=2, entities=[self.s_xneg]),
        )
        for material in materials:
            material.create_physical_group()

        boundary_groups = (
            BoundaryGroup(name="boundary_xneg", tag=10, dim=1, entities=[self.l_xneg]),
            BoundaryGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.l_xpos]),
            BoundaryGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.l_yneg0, self.l_yneg1]),
            BoundaryGroup(name="boundary_ypos", tag=13, dim=1, entities=[self.l_ypos0, self.l_ypos1]),
            BoundaryGroup(name="fault", tag=20, dim=1, entities=[self.l_fault]),
        )
        recursive = False
        if mark_vertices:
            boundary_groups += (BoundaryGroup(name="fault_end", tag=21, dim=0, entities=[self.p_fault_end]),
            )
            recursive = True
        for group in boundary_groups:
            group.create_physical_group(recursive=recursive)

    def generate_mesh(self):
        # right angle tris
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        if self.cell == "quad":
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)

        gmsh.model.mesh.generate(2)


class GenerateBox3D(GenerateApp):

    DOMAIN_X = DOMAIN_Y = DOMAIN_Z = 8.0e+3
    DX = 4.0e+3

    def __init__(self, cell="tet", mark_vertices=False, dump_mesh=False, gui=False):
        super().__init__(cell, mark_vertices, dump_mesh, gui)
        self.geometry = "box"

    def create_geometry(self):

        gmsh.initialize()
        gmsh.model.add("mesh")
        box = gmsh.model.occ.add_box(-0.5*self.DOMAIN_X, -0.5*self.DOMAIN_Y, -self.DOMAIN_Z,
                                    self.DOMAIN_X, self.DOMAIN_Y, self.DOMAIN_Z)

        disk = gmsh.model.occ.add_disk(0.0, 0.0, -0.5*self.DOMAIN_Z, self.DOMAIN_Y, self.DOMAIN_Z)
        gmsh.model.occ.rotate([(2,disk)], 0.0, 0.0, -0.5*self.DOMAIN_Z, 0.0, 1.0, 0.0, 0.5*numpy.pi)
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
        self.s_fault = gmsh.model.get_entities_in_bounding_box(-dx, bbox[1]-dx, bbox[2]-dx, +dx, bbox[4]+dx, bbox[5]+dx, dim=2)
        self.l_fault_end = gmsh.model.get_entities_in_bounding_box(-dx, bbox[1]-dx, bbox[2]-dx, +dx, bbox[1]+dx, bbox[5]+dx, dim=1)

        self.v_xneg = dimTags[0]
        self.v_xpos = dimTags[1]

    def mark(self, mark_vertices):
        materials = (
            MaterialGroup(tag=1, entities=[self.v_xneg[1]]),
            MaterialGroup(tag=2, entities=[self.v_xpos[1]]),
        )
        for material in materials:
            material.create_physical_group()

        boundary_groups = (
            BoundaryGroup(name="boundary_xneg", tag=10, dim=2, entities=[tag for dim, tag in self.s_xneg]),
            BoundaryGroup(name="boundary_xpos", tag=11, dim=2, entities=[tag for dim, tag in self.s_xpos]),
            BoundaryGroup(name="boundary_yneg", tag=12, dim=2, entities=[tag for dim, tag in self.s_yneg]),
            BoundaryGroup(name="boundary_ypos", tag=13, dim=2, entities=[tag for dim, tag in self.s_ypos]),
            BoundaryGroup(name="boundary_zneg", tag=14, dim=2, entities=[tag for dim, tag in self.s_zneg]),
            BoundaryGroup(name="boundary_zpos", tag=15, dim=2, entities=[tag for dim, tag in self.s_zpos]),
            BoundaryGroup(name="fault", tag=20, dim=2, entities=[tag for dim,tag in self.s_fault]),
        )
        recursive = False
        if mark_vertices:
            boundary_groups += (BoundaryGroup(name="fault_end", tag=21, dim=1, entities=[tag for dim,tag in self.l_fault_end]),
            )
            recursive = True
        for group in boundary_groups:
            group.create_physical_group(recursive=recursive)

    def generate_mesh(self):
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.DX)
        if self.cell == "hex":
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            
        gmsh.model.mesh.generate(3)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", action="store", dest="geometry", required=True, choices=["box-2d","layer-2d","box-3d"])
    parser.add_argument("--cell", action="store", dest="cell", required=True, choices=["tri","quad","hex","tet"])
    parser.add_argument("--dump-mesh", action="store_true", dest="dump_mesh")
    parser.add_argument("--gui", action="store_true", dest="gui")
    parser.add_argument("--mark-vertices", action="store_true", dest="mark_vertices")
    args = parser.parse_args()

    if args.geometry == "box-2d":
        app = GenerateBox2D
    elif args.geometry == "box-3d":
        app = GenerateBox3D
    app(args.cell, args.mark_vertices, args.dump_mesh, args.gui).main()
