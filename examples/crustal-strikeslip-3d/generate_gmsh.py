#!/usr/bin/env nemesis
"""Generate a tet mesh of 3 strike-slip faults using Gmsh, making use of the
Open Cascade geometry engine.

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
    #
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
    km = 1000.0
    DOMAIN_CENTER = (453700.0, 3947000.0)
    DOMAIN_X = 80*km
    DOMAIN_Y = 60*km
    DOMAIN_Z = 40*km
    FAULT_DEPTH = 15.0*km

    # Files with coordinates of fault traces in UTM zone 11
    FILENAME_MAINTRACE = "faulttrace_main_utm.txt"
    FILENAME_WESTTRACE = "faulttrace_west_utm.txt"
    FILENAME_STRANDETRACE = "faulttrace_east_utm.txt"

    # Discretization size on faults and bias (rate of cell size increase away from fault)
    DX_FAULT = 2.5*km
    DX_BIAS = 1.07

    def __init__(self):
        """Constructor.
        """
        super().__init__()

        # Set the cell choices available through command line options.
        # The default cell type `tri` and filename match the mesh used
        # in the PyLith parameter files.
        self.cell_choices = {
            "default": "tet",
            "choices": ["tet"],
            }
        self.filename = "mesh_tet.msh"

    def _create_points_from_file(self, filename):
        coordinates = numpy.loadtxt(filename)
        points = []
        for xy in coordinates:
            points.append(gmsh.model.occ.add_point(xy[0], xy[1], 0))
        return points

    def create_geometry(self):
        """Create geometry.

        We create the geometry in the UTM zone 11 coordinate system.

        We use the fault traces in UTM coordinates given in `faulttrace_{name}_utm.txt`.
        """
        # Create domain
        xSW = self.DOMAIN_CENTER[0] - 0.5*self.DOMAIN_X
        ySW = self.DOMAIN_CENTER[1] - 0.5*self.DOMAIN_Y
        zSW = -self.DOMAIN_Z
        self.v_domain = gmsh.model.occ.add_box(xSW, ySW, zSW, self.DOMAIN_X, self.DOMAIN_Y, self.DOMAIN_Z)

        # Create points for main fault
        points_fault_main = self._create_points_from_file(self.FILENAME_MAINTRACE)
        spline = gmsh.model.occ.add_spline(points_fault_main)
        dimTags = gmsh.model.occ.extrude([(1, spline)], 0, 0, -self.FAULT_DEPTH)
        self.s_fault_main = dimTags[1][1]

        # Create points for west fault strand
        points_fault_west = self._create_points_from_file(self.FILENAME_WESTTRACE)
        spline = gmsh.model.occ.add_spline(points_fault_west)
        dimTags = gmsh.model.occ.extrude([(1, spline)], 0, 0, -self.FAULT_DEPTH)
        self.s_fault_west = dimTags[1][1]

        # Create points for east fault strand
        points_fault_east = self._create_points_from_file(self.FILENAME_STRANDETRACE)
        spline = gmsh.model.occ.add_spline(points_fault_east)
        dimTags = gmsh.model.occ.extrude([(1, spline)], 0, 0, -self.FAULT_DEPTH)
        self.s_fault_east = dimTags[1][1]

        # Embed faults into domain so cells will align along fault
        gmsh.model.occ.fragment([(3, self.v_domain)], [(2, self.s_fault_main)], removeTool=True)
        gmsh.model.occ.fragment([(3, self.v_domain)], [(2, self.s_fault_west)], removeTool=True)
        gmsh.model.occ.fragment([(3, self.v_domain)], [(2, self.s_fault_east)], removeTool=True)
        gmsh.model.occ.synchronize()

        # Get surfaces of domain; identify order of surfaces using GUI
        _, surfaces = gmsh.model.get_adjacencies(3, self.v_domain)
        self.s_west = surfaces[0]
        self.s_south = surfaces[1]
        self.s_top = surfaces[2]
        self.s_north = surfaces[3]
        self.s_bottom = surfaces[4]
        self.s_east = surfaces[5]

        # Get fault surfaces and edges using the GUI
        self.s_fault_main_north = 17
        self.s_fault_main_south = 18
        self.s_fault_west = 16
        self.s_fault_east = 9

        self.fault_main_edges = [20, 21, 22, 23]
        self.fault_west_edges = [17, 18, 19]
        self.fault_east_edges = [19, 24, 25]

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Create a material for the domain.
        # The tag argument specifies the integer tag for the physical group.
        # The entities argument specifies the array of surfaces for the material.
        materials = (
            MaterialGroup(tag=1, entities=[self.v_domain]),
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries and the faults.
        # We use the `BoundaryGroup` data class defined in `gmsh_utils`.
        # The name and tag specify the name and tag assigned to the physical group.
        # The dimension and entities specify the geometric entities to include in the physical
        # group.
        face_groups = (
            BoundaryGroup(name="boundary_south", tag=10, dim=2, entities=[self.s_south]),
            BoundaryGroup(name="boundary_east", tag=11, dim=2, entities=[self.s_east]),
            BoundaryGroup(name="boundary_north", tag=12, dim=2, entities=[self.s_north]),
            BoundaryGroup(name="boundary_west", tag=13, dim=2, entities=[self.s_west]),
            BoundaryGroup(name="boundary_bottom", tag=14, dim=2, entities=[self.s_bottom]),
            BoundaryGroup(name="boundary_top", tag=15, dim=2, entities=[self.s_top]),

            BoundaryGroup(name="fault_main", tag=20, dim=2, entities=[self.s_fault_main_north, self.s_fault_main_south]),
            BoundaryGroup(name="fault_west", tag=21, dim=2, entities=[self.s_fault_west]),
            BoundaryGroup(name="fault_east", tag=22, dim=2, entities=[self.s_fault_east]),
            BoundaryGroup(name="fault_main_edges", tag=30, dim=1, entities=self.fault_main_edges),
            BoundaryGroup(name="fault_west_edges", tag=31, dim=1, entities=self.fault_west_edges),
            BoundaryGroup(name="fault_east_edges", tag=32, dim=1, entities=self.fault_east_edges),
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
        gmsh.model.mesh.field.setNumbers(field_distance, "SurfacesList", [self.s_fault_main_north, self.s_fault_main_south, self.s_fault_west, self.s_fault_east])

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

        if cell == "hex":
            # Generate a tri mesh and then recombine cells to form quadrilaterals.
            # We use the Frontal-Delaunay for Quads algorithm.
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.recombine()
        else:
            gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Laplace2D")


# If script is called from the command line, run the application.
if __name__ == "__main__":
    App().main()


# End of file

