#!/usr/bin/env nemesis
"""
We use the `gmsh_utils` module provided with PyLith. This module has helper functions
and classes for marking materials and boundaries compatible with PyLith. We also use the
`GenerateMesh` class from the module as a base class to our local `App` class. The
`GenerateMesh` class handles processing of command line options, initialization, and
finalizing the mesh.

Run `generate_gmsh.py --help` to see the command line options.

Run `generate_gmsh.py --write` to generate the mesh.
"""


# Import Gmsh Python interface
import gmsh

# Import the gmsh_utils Python module supplied with PyLith.
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Application used to generate the mesh using Gmsh.

    App uses `GenerateMesh` from `gmsh_utils` for common functionality that we avoid
    duplicating in each of our examples.
    """
    # Domain constants
    km = 1000.0
    LENGTH = 20*km
    DEPTH = 20*km

    # Magma conduit
    CONDUIT_LENGTH = 1*km
    CONDUIT_DEPTH = 6*km

    # Magma reservoir
    RES_CENTER = (LENGTH/2, -DEPTH+CONDUIT_DEPTH, 0)
    MAJOR_AXIS = 5*km
    MINOR_AXIS = 1*km

    DX_RES = 250
    DX_BIAS = 1.1

    def __init__(self):
        """Constructor.
        """
        super().__init__()

        # Set the cell choices available through command line options.
        # The default cell type `tri` and filename match the mesh used
        # in the PyLith parameter files.
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
            }
        self.filename = "mesh_tri.msh"

    def create_geometry(self):
        """Create geometry.
        """
        # Create domain surface
        domain_box = gmsh.model.occ.add_rectangle(0, -self.DEPTH, 0, self.LENGTH, self.DEPTH)
        gmsh.model.occ.add_plane_surface([domain_box])

        # Create chamber surface
        chamber = gmsh.model.occ.add_ellipse(self.RES_CENTER[0], self.RES_CENTER[1], self.RES_CENTER[2], self.MAJOR_AXIS, self.MINOR_AXIS)
        chamber_loop = gmsh.model.occ.add_curve_loop([chamber])
        s_chamber = gmsh.model.occ.add_plane_surface([chamber_loop])

        # Create conduit surface
        x0 = self.LENGTH/2 - self.CONDUIT_LENGTH/2
        y0 = -self.DEPTH
        s_conduit = gmsh.model.occ.add_rectangle(x0, y0, 0, self.CONDUIT_LENGTH, self.CONDUIT_DEPTH)

        # Join conduit and chamber into reservoir
        self.s_reservoir = gmsh.model.occ.fuse([(2, s_chamber)], [(2, s_conduit)])[0][0][1]

        # Embed reservoir in domain
        self.s_domain = gmsh.model.occ.fragment([(2, domain_box)], [(2, self.s_reservoir)])[0][1][1]

        gmsh.model.occ.remove_all_duplicates() # Duplicate points at intersection

        # Get boundaries
        _, domain_curves = gmsh.model.occ.get_curve_loops(self.s_domain)
        self.c_yneg_left = domain_curves[0][1]
        self.c_yneg_right = domain_curves[0][6]
        self.c_xpos = domain_curves[0][7]
        self.c_ypos = domain_curves[0][8]
        self.c_xneg = domain_curves[0][0]

        _, reservoir_curves = gmsh.model.occ.get_curve_loops(self.s_reservoir)
        self.c_flow = reservoir_curves[0][3]
        self.c_chamber_right = reservoir_curves[0][0]
        self.c_chamber_left = reservoir_curves[0][1]
        self.c_conduit_left = reservoir_curves[0][2]
        self.c_conduit_right = reservoir_curves[0][4]

        gmsh.model.occ.synchronize()

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
            MaterialGroup(tag=2, entities=[self.s_reservoir])
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries.
        # We use the `VertexGroup` data class defined in `gmsh_utils`.
        # The name and tag specify the name and tag assigned to the physical group.
        # The dimension and entities specify the geometric entities to include in the physical
        # group.
        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=20, dim=1, entities=[self.c_xneg]),
            VertexGroup(name="boundary_xpos", tag=21, dim=1, entities=[self.c_xpos]),
            VertexGroup(name="boundary_yneg", tag=22, dim=1, entities=[self.c_yneg_left, self.c_flow, self.c_yneg_right]),
            VertexGroup(name="boundary_ypos", tag=23, dim=1, entities=[self.c_ypos]),
            VertexGroup(name="boundary_flow", tag=24, dim=1, entities=[self.c_flow])
        )
        for group in vertex_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Set discretization size with geometric progression from distance to the reservoir boundary.

        # We turn off the default sizing methods.
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # First, we setup a field `field_distance` with the distance from the reservoir boundary.
        distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(distance, "CurvesList", [self.c_chamber_right, self.c_chamber_left, self.c_conduit_left, self.c_flow, self.c_conduit_right])

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the boundary, the distance from
        # the boundary (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression`
        # for creating the string with the mathematical function.
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(distance, min_dx=self.DX_RES, bias=self.DX_BIAS)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        # Finally, we use the field `field_size` for the cell size of the mesh.
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        if cell == "quad":
            # Generate a tri mesh and then recombine cells to form quadrilaterals.
            # We use the Frontal-Delaunay for Quads algorithm.
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.recombine()
            gmsh.model.mesh.generate(2)
        else:
            gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")


# If script is called from the command line, run the application.
if __name__ == "__main__":
    App().main()


# End of file

