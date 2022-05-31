#!/usr/bin/env nemesis
"""Generate a tri or quad mesh of a reverse fault with a splay fault using Gmsh, making
use of the built-in geometry engine.

We use the `gmsh_utils` module provided with PyLith. This module has helper functions
and classes for marking materials and boundaries compatible with PyLith. We also use the
`GenerateMesh` class from the module as a base class to our local `App` class. The
`GenerateMesh` class handles processing of command line options, initialization, and
finalizing the mesh.

Run `generate_gmsh.py --help` to see the command line options.

Run `generate_gmsh.py --write` to generate the mesh.
"""

# Import the math module for trigonometric functions
import math

# Import Gmsh Python interface
import gmsh

# Import the gmsh_utils Python module supplied with PyLith.
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Application used to generate the mesh using Gmsh.

    App uses `GenerateMesh` from `gmsh_utils` for common functionality that we avoid
    duplicating in each of our examples.

    Domain is 200km by 100km.

    -100.0 km <= x <= 100.0 km
    -100.0 km <= y <= 0.0 km

    We create a reverse fault and a splay fault.
    """
    # Define some constants that determine the geometry of the domain.
    DOMAIN_X = 200.0e+3
    DOMAIN_Y = 100.0e+3
    FAULT_WIDTH = 60.0e+3
    FAULT_DIP = 30.0
    SPLAY_DIP = 45.0
    SPLAY_OFFSET = 15.0e+3

    def __init__(self):
        """Constructor.
        """
        # Set the cell choices available through command line options
        # with the default cell type `tri` matching the PyLith parameter files.
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
            }
        self.filename = "mesh_tri.msh"

    def create_geometry(self):
        """Create geometry.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Set local variables for domain size and corner of the domain.
        lx = self.DOMAIN_X
        ly = self.DOMAIN_Y
        x1 = -0.5 * lx
        y1 = -ly

        # Create points for domain perimeter.
        p1 = gmsh.model.geo.add_point(x1, y1, 0.0)
        p2 = gmsh.model.geo.add_point(x1+lx, y1, 0.0)
        p3 = gmsh.model.geo.add_point(x1+lx, y1+ly, 0.0)
        p4 = gmsh.model.geo.add_point(x1, y1+ly, 0.0)

        # Create points for main fault
        w = self.FAULT_WIDTH
        fault_dip = self.FAULT_DIP / 180.0 * math.pi
        x6 = -0.5 * w * math.cos(fault_dip)
        y6 = -w * math.sin(fault_dip)
        x7 = x6+w*math.cos(fault_dip)
        y5 = -(x7-x1)*math.tan(fault_dip)
        p5 = gmsh.model.geo.add_point(x1, y5, 0.0)
        p6 = gmsh.model.geo.add_point(x6, y6, 0.0)
        p7 = gmsh.model.geo.add_point(x7, 0.0, 0.0)
        self.p_fault_end = p6

        # Create points for splay fault
        splay_dip = self.SPLAY_DIP / 180.0 * math.pi
        w = self.SPLAY_OFFSET * math.sin(fault_dip) / (math.sin(splay_dip)*math.cos(fault_dip) - math.sin(fault_dip)*math.cos(splay_dip))
        x9 = x7 - self.SPLAY_OFFSET
        x8 = x9 - w*math.cos(splay_dip)
        y8 = -w * math.sin(splay_dip)
        p8 = gmsh.model.geo.add_point(x8, y8, 0.0)
        p9 = gmsh.model.geo.add_point(x9, 0.0, 0.0)
        self.p_splay_end = p8

        # Create curves. We store the curve tag as a data member
        # so that we can refer to them later.
        # We traverse the curves in a counter clock-wise direction.
        # If the curve is in the opporite direction, we use the negative tag.
        self.c_yneg = gmsh.model.geo.add_line(p1, p2)
        self.c_xpos = gmsh.model.geo.add_line(p2, p3)
        self.c_ypos_fw = gmsh.model.geo.add_line(p3, p7)
        self.c_ypos_w = gmsh.model.geo.add_line(p7, p9)
        self.c_ypos_hw = gmsh.model.geo.add_line(p9, p4)
        self.c_xneg_hw = gmsh.model.geo.add_line(p4, p5)
        self.c_xneg_fw = gmsh.model.geo.add_line(p5, p1)

        self.c_fault_ext = gmsh.model.geo.add_line(p5, p6)
        self.c_fault_l = gmsh.model.geo.add_line(p6, p8)
        self.c_fault_u = gmsh.model.geo.add_line(p8, p7)

        self.c_splay = gmsh.model.geo.add_line(p8, p9)

        c0 = gmsh.model.geo.add_curve_loop([self.c_yneg, self.c_xpos, self.c_ypos_fw, -self.c_fault_u, -self.c_fault_l, -self.c_fault_ext, self.c_xneg_fw])
        self.s_slab = gmsh.model.geo.add_plane_surface([c0])
        c0 = gmsh.model.geo.add_curve_loop([self.c_fault_ext, self.c_fault_l, self.c_splay, self.c_ypos_hw, self.c_xneg_hw])
        self.s_plate = gmsh.model.geo.add_plane_surface([c0])
        c0 = gmsh.model.geo.add_curve_loop([self.c_fault_u, self.c_ypos_w, -self.c_splay])
        self.s_wedge = gmsh.model.geo.add_plane_surface([c0])

        gmsh.model.geo.synchronize()


    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Create three materials (slab, plate, and wedge).
        # We use the `MaterialGroup` data class defined in `gmsh_utils.`
        # The tag argument specifies the integer tag for the physical group.
        # The entities argument specifies the array of surfaces for the material.
        materials = (
            MaterialGroup(tag=1, entities=[self.s_slab]),
            MaterialGroup(tag=2, entities=[self.s_plate]),
            MaterialGroup(tag=3, entities=[self.s_wedge]),
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries and the faults.
        # We use the `VertexGroup` data class defined in `gmsh_utils`.
        # The name and tag specify the name and tag assigned to the physical group.
        # The dimension and entities specify the geometric entities to include in the physical
        # group.
        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=10, dim=1, entities=[self.c_xneg_hw, self.c_xneg_fw]),
            VertexGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.c_xpos]),
            VertexGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.c_yneg]),
            VertexGroup(name="boundary_ypos", tag=13, dim=1, entities=[self.c_ypos_fw, self.c_ypos_w, self.c_ypos_hw]),
            VertexGroup(name="fault", tag=20, dim=1, entities=[self.c_fault_l, self.c_fault_u]),
            VertexGroup(name="fault_end", tag=21, dim=0, entities=[self.p_fault_end]),
            VertexGroup(name="splay", tag=22, dim=1, entities=[self.c_splay]),
            VertexGroup(name="splay_end", tag=23, dim=0, entities=[self.p_splay_end]),
        )
        for group in vertex_groups:
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
        fault_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(fault_distance, "CurvesList", [self.c_fault_l, self.c_fault_u])

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the fault, the distance from
        # the fault (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression` 
        # for creating the string with the mathematical function.
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(fault_distance, min_dx=2.0e+3, bias=1.05)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        ## Finally, we use the field `field_size` for the cell size of the mesh.
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
