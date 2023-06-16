"""
Utilities for generating PyLith compatible meshes using Gmsh.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass

import gmsh

@dataclass
class VertexGroup:
    name: str
    tag: int
    dim: int
    entities: list

    def create_physical_group(self, recursive=True):
        create_group(self.name, self.tag, self.dim, self.entities, recursive)

@dataclass
class MaterialGroup:
    tag: int
    entities: list

    def create_physical_group(self):
        create_material(self.tag, self.entities)


def create_group(name, tag, dim, entities, recursive=True, exclude=None):
    gmsh.model.add_physical_group(dim, entities, tag)
    gmsh.model.set_physical_name(dim, tag, name)
    entities_lowerdim = []
    for entity in entities:
        entities_up, entities_down = gmsh.model.get_adjacencies(dim, entity)
        entities_lowerdim += [e for e in entities_down]
    if recursive and dim >= 1:
        create_group(name, tag, dim-1, entities_lowerdim)


def get_physical_group(name):
    """Get all physical groups matching name.
    """
    groups = {}
    dimTags = gmsh.model.get_physical_groups()
    for dim, tag in dimTags:
        if gmsh.model.get_physical_name(dim, tag) == name:
            groups[dim] = tag
    if len(groups) == 0:
        raise ValueError(f"Could not find physical group '{name}'.")
    return groups

def group_exclude(group_name, exclude_name, new_name, new_tag):
    """Remove entities from physical group `exclude_name` from physical group `group_name`.
    """
    groups = get_physical_group(group_name)
    groups_exclude = get_physical_group(exclude_name)

    for dim, tag in groups.items():
        entities = set(gmsh.model.get_entities_for_physical_group(dim, tag))
        if dim in groups_exclude:
            entities_exclude = gmsh.model.get_entities_for_physical_group(dim, groups_exclude[dim])
            entities = entities.difference(set(entities_exclude))
        tags = [tag for tag in entities]
        gmsh.model.add_physical_group(dim, tags, new_tag)
        gmsh.model.set_physical_name(dim, new_tag, new_name)


def create_material(tag, entities):
    if tag <= 0:
        raise ValueError(f"ERROR: Attempting to use non-positive material tag '{tag}'. Tags for physical groups must be positive.")
    dim = gmsh.model.get_dimension()
    name = gmsh.model.get_physical_name(dim, tag)
    if name:
        raise ValueError(f"ERROR: Attempting to use material tag '{tag}' that is already in use for material '{name}'.")
    gmsh.model.addPhysicalGroup(dim, entities, tag)
    gmsh.model.setPhysicalName(dim, tag, f"material-id:{tag}")


class GenerateMesh(ABC):

    @staticmethod
    def get_math_progression(field_distance, min_dx, bias):
        """Generate the Gmsh MathEval string corresponding to the cell size as a function
        of distance, starting cell size, and bias factor.

        The expression is min_dx * bias**n, where n is the number of cells from the fault.
        n = log(1+distance/min_dx*(bias-1))/log(bias)

        In finding the expression for `n`, we make use that the sum of a geometric series with n
        terms Sn = min_dx * (1 + bias + bias**2 + ... + bias**n) = min_dx * (bias**n - 1)/(bias - 1).
        """
        return f"{min_dx}*{bias}^(Log(1.0+F{field_distance}/{min_dx}*({bias}-1.0))/Log({bias}))"

    def __init__(self):
        self.cell_choices = {"required": False}
        self.filename = "mesh.msh"

    def main(self):
        """
        Main entry point for meshing application.
        """
        args = self._parse_command_line()

        self.initialize(args.name)
        if args.geometry:
            self.create_geometry()
        if args.mark:
            self.mark()
        if args.generate:
            self.generate_mesh(args.cell)
        if args.write:
            self.write(args.filename, args.binary)
        self.finalize(args.gui)

    def initialize(self, name: str):
        """Initialize Gmsh.

        :param name: Name for mesh.
        """
        gmsh.initialize()
        gmsh.model.add(name)

    def finalize(self, gui=False):
        """Finalize Gmsh.

        :param gui: Show GUI if True.
        """
        if gui:
            gmsh.fltk.run()
        gmsh.finalize()
    
    def write(self, filename: str, binary=True):
        """Write mesh to file.

        :param filename: Name of file.
        :param binary: Write as binary file if True.
        """
        binary_flag = 1 if binary else 0
        gmsh.option.setNumber("Mesh.Binary", binary_flag)
        gmsh.write(filename)

    @abstractmethod
    def create_geometry(self):
        """Create geometry.
        """
        pass

    @abstractmethod
    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.
        """
        pass

    @abstractmethod
    def generate_mesh(self, cell):
        """Generate the mesh. Should also include optimizing the mesh quality.
        """
        pass

    def _parse_command_line(self):
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--geometry", action="store_true", dest="geometry", help="Create geometry.")
        parser.add_argument("--mark", action="store_true", dest="mark", help="Mark entities (materials, boundary conditions).")
        parser.add_argument("--generate", action="store_true", dest="generate", help="Generate mesh.")
        parser.add_argument("--write", action="store_true", dest="write", help="Write mesh to file.")

        parser.add_argument("--name", action="store", dest="name", help="Name of mesh in Gmsh.", default="mesh")
        parser.add_argument("--filename", action="store", dest="filename", help="Name of output mesh file.", default=self.filename)
        parser.add_argument("--ascii", action="store_false", dest="binary", help="Write mesh to ASCII file (detault is binary).", default=True)

        parser.add_argument("--cell", action="store", dest="cell", **self.cell_choices)

        parser.add_argument("--gui", action="store_true", dest="gui", help="Show GUI after running steps.")

        GenerateMesh._add_arguments(parser)
        args = parser.parse_args()
        if args.write:
            args.generate = True
        if args.generate:
            args.mark = True
        if args.mark:
            args.geometry = True
        return args

    @staticmethod
    def _add_arguments(parser):
        pass


# End of file

