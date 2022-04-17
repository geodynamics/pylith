#!/usr/bin/env nemesis
"""
Physical groups have a name, dimension, and physical tag. We will make a label for each physical group.
Entities have a dimension, bounding box, and physical tag
Nodes have a node tag, coordinates, and an entity. We need to lookup the entity to get the physical tag for a node.
Elements have an element tag, dimension, element type, list of node tags, and an entity. We need to lookup the entity to get the physical tag for a node.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
import gmsh
import numpy

def create_group(name, tag, dim, entities, recursive=True):
    gmsh.model.add_physical_group(dim, entities, tag)
    gmsh.model.set_physical_name(dim, tag, name)
    entities_lowerdim = []
    for entity in entities:
        entities_up, entities_down = gmsh.model.get_adjacencies(dim, entity)
        entities_lowerdim += [e for e in entities_down]
    if recursive and dim >= 1:
        create_group(name, tag, dim-1, entities_lowerdim)


def create_material(tag, entities):
    if tag <= 0:
        raise ValueError(f"ERROR: Attempting to use non-positive material tag '{tag}'. Tags for physical groups must be positive.")
    dim = gmsh.model.get_dimension()
    name = gmsh.model.get_physical_name(dim, tag)
    if name:
        raise ValueError(f"ERROR: Attempting to use material tag '{tag}' that is already in use for material '{name}'.")
    gmsh.model.addPhysicalGroup(dim, entities, tag)
    gmsh.model.setPhysicalName(dim, tag, f"material-id:{tag}")


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
        entities = gmsh.model.mesh.get_elements(dim, tag)
        cell_indices.append(entities[1][0])
        min_tag = min(numpy.min(entities[1][0]), min_tag)
    for indices in cell_indices:
        indices -= min_tag
    material_ids = numpy.zeros(cells.shape[0], dtype=numpy.uint64)
    for (dim, tag), indices in zip(physical_groups, cell_indices):
        material_ids[indices] = tag
    return material_ids


def get_vertex_groups():
    dim = gmsh.model.get_dimension()-1
    physical_groups = gmsh.model.get_physical_groups(dim)
    groups = {}
    for dim, tag in physical_groups:
        name = gmsh.model.get_physical_name(dim, tag)
        node_tags, node_coords = gmsh.model.mesh.get_nodes_for_physical_group(dim, tag)
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
        self.fout.write(f"static const PylithScalar vertices[{nvertices}*{space_dim}] = " + "{\n")
        for vertex in vertices:
            line = ", ".join([f"{x:+16.8e}" for x in vertex]) + ",\n"
            self.fout.write(line)
        self.fout.write("};\n")
        self.fout.write("_data->vertices = const_cast<PylithScalar*>(vertices);\n\n")

    def write_cells(self, cells):
        ncells, ncorners = cells.shape
        print(f"static const PylithInt cells[{ncells}*{ncorners}] = " + "{")
        for cell in cells[:,0:ncorners]:
            line = ", ".join([f"{int(node-1):d}" for node in cell]) + ",\n"
            self.fout.write(line)
        self.fout.write("};\n")
        self.fout.write("_data->cells = const_cast<PylithInt*>(cells);\n\n")

    def write_material_ids(self, material_ids):
        ncells = len(material_ids)
        self.fout.write(f"static const PylithInt materialIds[{ncells}] = " + "{\n")
        line = [f"{mid}," for mid in material_ids]
        self.fout.write("\n".join(line))
        self.fout.write("\n};\n")
        self.fout.write("_data->materialIds = const_cast<PylithInt*>(materialIds);\n\n")

    def write_vertex_groups(self, groups):
        ngroups = len(groups)
        self.fout.write(f"_data->numGroups = {ngroups};\n")
        gsizes = [len(group) for group in groups.values()]
        gsizes_str = [f"{len(group)}" for group in groups.values()]
        self.fout.write(f"static const PylithInt groupSizes[{ngroups}] = ")
        self.fout.write("{" + ", ".join(gsizes_str) + "};\n")
        self.fout.write("_data->groupSizes = const_cast<PylithInt*>(groupSizes);\n")
        self.fout.write(f"static const PylithInt groups[{sum(gsizes)}] = " + "{\n")
        for group in groups.values():
            line = ", ".join([f"{int(vertex-1):d}" for vertex in group]) + ",\n"
            self.fout.write(line)
        self.fout.write("};\n")
        self.fout.write("_data->groups = const_cast<PylithInt*>(groups);\n")
        self.fout.write(f"static const char* groupNames[{ngroups}] = " + "{\n")
        for name in groups.keys():
            self.fout.write(f'"{name}",\n')
        self.fout.write("};\n")
        self.fout.write("_data->groupNames = const_cast<char**>(groupNames);\n")

def write_mesh(writer=WriterCxx()):
    writer.write_vertices(get_vertices())
    cells = get_cells()
    writer.write_cells(cells)
    writer.write_material_ids(get_material_ids(cells))
    del cells
    writer.write_vertex_groups(get_vertex_groups())

@dataclass
class VertexGroup:
    name: str
    tag: int
    entities: list

@dataclass
class MaterialGroup:
    tag: int
    entities: list

class GenerateApp(ABC):

    def __init__(self, cell, dump_mesh=False, gui=False):
        self.cell = cell
        self.dump_mesh = dump_mesh
        self.gui = gui

    def main(self):
        self.create_geometry()
        self.mark()
        self.generate_mesh()
        self.write()
        if self.dump_mesh:
            self.print_mesh()
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
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(f"{self.cell}_ascii.msh")

        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(f"{self.cell}_binary.msh")

    def print_mesh(self):
        write_mesh(writer=WriterCxx())

    def finalize(self):
        if self.gui:
            gmsh.fltk.run()
        gmsh.finalize()


class Generate2D(GenerateApp):
    """
    v5----v4----v3
    |      |     |
    |      |     |
    |      |     |
    v0----v1----v2
    """

    DOMAIN_X = DOMAIN_Y = 8.0e+3
    DX = 4.0e+3

    def __init__(self, cell="tri", dump_mesh=False, gui=False):
        super().__init__(cell, dump_mesh, gui)

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

        c0 = gmsh.model.geo.addCurveLoop([self.l_yneg0, self.l_fault, self.l_ypos1, self.l_xneg])
        self.s_xneg = gmsh.model.geo.addPlaneSurface([c0])
        c1 = gmsh.model.geo.addCurveLoop([self.l_yneg1, self.l_xpos, self.l_ypos0, -self.l_fault])
        self.s_xpos = gmsh.model.geo.addPlaneSurface([c1])

        gmsh.model.geo.synchronize()

    def mark(self):
        materials = (
            MaterialGroup(tag=1, entities=[self.s_xneg]),
            MaterialGroup(tag=2, entities=[self.s_xpos]),
        )
        for material in materials:
            create_material(material.tag, material.entities)

        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=10, entities=[self.l_xneg]),
            VertexGroup(name="boundary_xpos", tag=11, entities=[self.l_xpos]),
            VertexGroup(name="boundary_yneg", tag=12, entities=[self.l_yneg0, self.l_yneg1]),
            VertexGroup(name="boundary_ypos", tag=13, entities=[self.l_ypos0, self.l_ypos1]),
            VertexGroup(name="fault", tag=20, entities=[self.l_fault]),
        )
        for group in vertex_groups:
            create_group(group.name, group.tag, gmsh.model.get_dimension()-1, group.entities)

    def generate_mesh(self):
        # right angle tris
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        if self.cell == "quad":
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)

        gmsh.model.mesh.generate(2)


class Generate3D(GenerateApp):

    DOMAIN_X = DOMAIN_Y = DOMAIN_Z = 8.0e+3
    DX = 4.0e+3

    def __init__(self, cell="tet", dump_mesh=False, gui=False):
        super().__init__(cell, dump_mesh, gui)

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

        self.v_xneg = dimTags[0]
        self.v_xpos = dimTags[1]

    def mark(self):
        materials = (
            MaterialGroup(tag=1, entities=[self.v_xneg[1]]),
            MaterialGroup(tag=2, entities=[self.v_xpos[1]]),
        )
        for material in materials:
            create_material(material.tag, material.entities)

        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=10, entities=[tag for dim, tag in self.s_xneg]),
            VertexGroup(name="boundary_xpos", tag=11, entities=[tag for dim, tag in self.s_xpos]),
            VertexGroup(name="boundary_yneg", tag=12, entities=[tag for dim, tag in self.s_yneg]),
            VertexGroup(name="boundary_ypos", tag=13, entities=[tag for dim, tag in self.s_ypos]),
            VertexGroup(name="boundary_zneg", tag=14, entities=[tag for dim, tag in self.s_zneg]),
            VertexGroup(name="boundary_zpos", tag=15, entities=[tag for dim, tag in self.s_zpos]),
            VertexGroup(name="fault", tag=20, entities=[tag for dim,tag in self.s_fault]),
        )
        for group in vertex_groups:
            create_group(group.name, group.tag, gmsh.model.get_dimension()-1, group.entities)

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
    parser.add_argument("--cell", action="store", dest="cell", required=True, choices=["tri","quad","hex","tet"])
    parser.add_argument("--dump-mesh", action="store_true", dest="dump_mesh")
    parser.add_argument("--gui", action="store_true", dest="gui")
    args = parser.parse_args()

    if args.cell in ["tri", "quad"]:
        app = Generate2D
    elif args.cell in ["tet", "hex"]:
        app = Generate3D
    app(args.cell, args.dump_mesh, args.gui).main()
