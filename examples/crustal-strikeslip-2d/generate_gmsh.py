#!/usr/bin/env nemesis
"""Generate a tri of 3 strike-slip faults using Gmsh, making use of the
built-in geometry engine.

We use the `gmsh_utils` module provided with PyLith. This module has helper functions
and classes for marking materials and boundaries compatible with PyLith. We also use the
`GenerateMesh` class from the module as a base class to our local `App` class. The
`GenerateMesh` class handles processing of command line options, initialization, and
finalizing the mesh.

Run `generate_gmsh.py --help` to see the command line options.

Run `generate_gmsh.py --write` to generate the mesh.
"""

import numpy

# Import Gmsh Python interface
import gmsh

# Import the gmsh_utils Python module supplied with PyLith.
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Application used to generate the mesh using Gmsh.

    App uses `GenerateMesh` from `gmsh_utils` for common functionality that we avoid
    duplicating in each of our examples.
    """
    #    pNW---------------------------------pNE
    #    |                                   |
    #    |                                   |
    #    |         main                      |
    #    |          \          east          |
    #    |           \        /              |
    #    |            \      /               |
    #    |             \    /                |
    #    |              pI                   |
    #    |             /   \                 |
    #    |            /     \                |
    #    |          /        \               |
    #    |        /                          |
    #    |      /                            |
    #    |     west                          |
    #    |                                   |
    #    |                                   |
    #    pSW---------------------------------pSE
    #
    #    The fault traces intersect at point pI.

    # Locations defining domain
    DOMAIN_W = 410000.00
    DOMAIN_E = 490000.00
    DOMAIN_S = 3910000.00
    DOMAIN_N = 3990000.00

    # Files with coordinates of fault traces in UTM zone 11
    FILENAME_MAINTRACE = "faulttrace_main_utm.txt"
    FILENAME_WESTTRACE = "faulttrace_west_utm.txt"
    FILENAME_STRANDETRACE = "faulttrace_east_utm.txt"

    # Discretization size on faults and bias (rate of cell size increase away from fault)
    DX_FAULT = 1.0e+3
    DX_BIAS = 1.07

    def __init__(self):
        """Constructor.
        """
        super().__init__()

        # Set the cell choices available through command line options.
        # The default cell type `tri` and filename match the mesh used
        # in the PyLith parameter files.
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri", "quad"],
            }
        self.filename = "mesh_tri.msh"

    def _create_points_from_file(self, filename):
        coordinates = numpy.loadtxt(filename)
        points = []
        for xy in coordinates:
            points.append(gmsh.model.geo.add_point(xy[0], xy[1], 0))
        return points

    def create_geometry(self):
        """Create geometry.

        We create the geometry in the UTM zone 11 coordinate system.

        We use the fault traces in UTM coordinates given in `faulttrace_{name}_utm.txt`.
        """
        # Create points for domain
        pSW = gmsh.model.geo.add_point(self.DOMAIN_W, self.DOMAIN_S, 0)
        pSE = gmsh.model.geo.add_point(self.DOMAIN_E, self.DOMAIN_S, 0)
        pNE = gmsh.model.geo.add_point(self.DOMAIN_E, self.DOMAIN_N, 0)
        pNW = gmsh.model.geo.add_point(self.DOMAIN_W, self.DOMAIN_N, 0)

        # Create curves for domain
        self.c_south = gmsh.model.geo.add_line(pSW, pSE)
        self.c_east = gmsh.model.geo.add_line(pSE, pNE)
        self.c_north = gmsh.model.geo.add_line(pNE, pNW)
        self.c_west = gmsh.model.geo.add_line(pNW, pSW)

        # Create curve loop and surface for domain
        loop = gmsh.model.geo.add_curve_loop([self.c_south, self.c_east, self.c_north, self.c_west])
        self.s_domain = gmsh.model.geo.add_plane_surface([loop])


        # Create points for main fault
        points = self._create_points_from_file(self.FILENAME_MAINTRACE)
        p_fault_intersection = points[3] # Intersection of fault strands and main fault
        spline = gmsh.model.geo.add_spline(points)
        self.curves_fault_main = gmsh.model.geo.split_curve(spline, [p_fault_intersection])
        self.fault_main_ends = [points[0], points[-1]] # Ends of fault

        # Create points for west fault strand
        points = self._create_points_from_file(self.FILENAME_WESTTRACE)
        self.c_fault_west = gmsh.model.geo.add_spline(points)

        # Create points for east fault strand
        points = self._create_points_from_file(self.FILENAME_STRANDETRACE)
        self.c_fault_east = gmsh.model.geo.add_spline(points)


        # Embed faults into domain so cells will align along fault
        gmsh.model.geo.remove_all_duplicates() # Duplicate points at intersection
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.embed(1, self.curves_fault_main, 2, self.s_domain)
        gmsh.model.mesh.embed(1, [self.c_fault_west], 2, self.s_domain)
        gmsh.model.mesh.embed(1, [self.c_fault_east], 2, self.s_domain)
        gmsh.model.geo.synchronize()

        # Get points at ends of branches after removing duplicate points
        _, self.fault_west_ends = gmsh.model.get_adjacencies(dim=1, tag=self.c_fault_west)
        _, self.fault_east_ends = gmsh.model.get_adjacencies(dim=1, tag=self.c_fault_east)


    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Create a material for the domain.
        # The tag argument specifies the integer tag for the physical group.
        # The entities argument specifies the array of surfaces for the material.
        materials = (
            MaterialGroup(tag=1, entities=[self.s_domain]),
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries and the faults.
        # We use the `BoundaryGroup` data class defined in `gmsh_utils`.
        # The name and tag specify the name and tag assigned to the physical group.
        # The dimension and entities specify the geometric entities to include in the physical
        # group.
        face_groups = (
            BoundaryGroup(name="boundary_south", tag=10, dim=1, entities=[self.c_south]),
            BoundaryGroup(name="boundary_east", tag=11, dim=1, entities=[self.c_east]),
            BoundaryGroup(name="boundary_north", tag=12, dim=1, entities=[self.c_north]),
            BoundaryGroup(name="boundary_west", tag=13, dim=1, entities=[self.c_west]),

            BoundaryGroup(name="fault_main", tag=20, dim=1, entities=self.curves_fault_main),
            BoundaryGroup(name="fault_west", tag=21, dim=1, entities=[self.c_fault_west]),
            BoundaryGroup(name="fault_east", tag=22, dim=1, entities=[self.c_fault_east]),
            BoundaryGroup(name="fault_main_ends", tag=30, dim=0, entities=self.fault_main_ends),
            BoundaryGroup(name="fault_west_ends", tag=31, dim=0, entities=self.fault_west_ends),
            BoundaryGroup(name="fault_east_ends", tag=32, dim=0, entities=self.fault_east_ends),
        )
        for group in face_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Set discretization size with geometric progression from distance to the fault.
        
        # We turn off the default sizing methods.
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # First, we setup a field `field_distance` with the distance from the fault.
        field_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_distance, "CurvesList", self.curves_fault_main.tolist() + [self.c_fault_west, self.c_fault_east])

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the fault, the distance from
        # the fault (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression` 
        # for creating the string with the mathematical function.
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(field_distance, min_dx=self.DX_FAULT, bias=self.DX_BIAS)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        # Finally, we use the field `field_size` for the cell size of the mesh.
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        if cell == "quad":
            # Generate a tri mesh and then recombine cells to form quadrilaterals.
            # We use the Frontal-Delaunay for Quads algorithm.
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.recombine()
        else:
            gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")


# If script is called from the command line, run the application.
if __name__ == "__main__":
    App().main()


# End of file

