#!/usr/bin/env nemesis
"""Generate a tri or quad mesh of a strike-slip fault using Gmsh, making use of the
built-in geometry engine.

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
import numpy as np
import pandas as pd
# Import the gmsh_utils Python module supplied with PyLith.
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Application used to generate the mesh using Gmsh.

    App uses `GenerateMesh` from `gmsh_utils` for common functionality that we avoid
    duplicating in each of our examples.

    Domain is 100km by 150km.
    -50.0 km <= x <= 50.0 km
    -75.0 km <= y <= 75.0 km

    The fault surface runs along the y-axis through the entire domain.

    p4-----p6-----p3
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    p1-----p5-----p2
    """
    #xy = pd.read_csv('xy.txt', sep = ' ', names = ['x','y'])
    

    def __init__(self):
        """Constructor.
        """
        # Set the cell choices available through command line options.
        # The default cell type `tri` and filename match the mesh used
        # in the PyLith parameter files.
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
            }
        self.filename = "mesh_tri_2.msh"

    def create_geometry(self):
        """Create geometry.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Set local variables for domain size and corner of the domain.
       

        # Create points.
        #p1 = gmsh.model.geo.add_point(xy, y1, 0.0)
        xy = pd.read_csv('xy.txt', sep = ' ', names = ['x','y'])
        points = []
        for i in range(len(xy.x)):
            points.append(gmsh.model.geo.add_point(xy.x[i],xy.y[i],0))
        self.p_inter = gmsh.model.geo.add_point(xy.x[7],xy.y[7],0)

        # Create curves. We store the curve tag as a data member
        # so that we can refer to them later.
        self.c_yneg = gmsh.model.geo.add_line(points[0], points[1])
        self.c_xpos = gmsh.model.geo.add_line(points[1], points[2])
        self.c_ypos = gmsh.model.geo.add_line(points[3], points[2])
        self.c_xneg = gmsh.model.geo.add_line(points[0], points[3])

        self.c_link1 = gmsh.model.geo.add_line(points[0], points[11])
        self.c_link2 = gmsh.model.geo.add_line(points[1], points[10])
        self.c_link3 = gmsh.model.geo.add_line(points[17], points[2])
        self.c_link4 = gmsh.model.geo.add_line(points[4], points[3])

        self.c_fault1sge1 = gmsh.model.geo.add_line(points[10], points[9])
        self.c_fault1sge2 = gmsh.model.geo.add_line(points[9], points[8])
        self.c_fault1sge3 = gmsh.model.geo.add_line(points[8], points[7])
        self.c_fault1sge4 = gmsh.model.geo.add_line(points[7], points[6])
        self.c_fault1sge5 = gmsh.model.geo.add_line(points[6], points[5])
        self.c_fault1sge6 = gmsh.model.geo.add_line(points[5], points[4])
        
        self.c_fault2sge1 = gmsh.model.geo.add_line(points[11], points[12])
        self.c_fault2sge2 = gmsh.model.geo.add_line(points[12], points[13])
        self.c_fault2sge3 = gmsh.model.geo.add_line(points[13], points[14])
        self.c_fault2sge4 = gmsh.model.geo.add_line(points[14], points[7])
        self.c_fault2sge5 = gmsh.model.geo.add_line(points[7], points[15])
        self.c_fault2sge6 = gmsh.model.geo.add_line(points[15], points[16])
        self.c_fault2sge7 = gmsh.model.geo.add_line(points[16], points[17])
      


        # Create curve loops and surfaces from the curves.
        # We traverse the curves in a counter clock-wise direction.
        # If the curve is in the opporite direction, we use the negative tag.
        c0 = gmsh.model.geo.add_curve_loop([self.c_link1, self.c_fault2sge1, self.c_fault2sge2,self.c_fault2sge3,self.c_fault2sge4,self.c_fault1sge4,self.c_fault1sge5,self.c_fault1sge6,self.c_link4,-self.c_xneg])
        self.s_1 = gmsh.model.geo.add_plane_surface([c0])
        c1 = gmsh.model.geo.add_curve_loop([self.c_yneg,self.c_link2,self.c_fault1sge1,self.c_fault1sge2,self.c_fault1sge3,-self.c_fault2sge4,-self.c_fault2sge3,-self.c_fault2sge2,-self.c_fault2sge1,-self.c_link1])
        self.s_2 = gmsh.model.geo.add_plane_surface([c1])
        c2 = gmsh.model.geo.add_curve_loop([self.c_xpos,-self.c_link3,-self.c_fault2sge7,-self.c_fault2sge6,-self.c_fault2sge5,-self.c_fault1sge3,-self.c_fault1sge2,-self.c_fault1sge1,-self.c_link2])
        self.s_3 = gmsh.model.geo.add_plane_surface([c2])
        c3 = gmsh.model.geo.add_curve_loop([-self.c_ypos,-self.c_link4,-self.c_fault1sge6,-self.c_fault1sge5,-self.c_fault1sge4,self.c_fault2sge5,self.c_fault2sge6,self.c_fault2sge7,self.c_link3])
        self.s_4 = gmsh.model.geo.add_plane_surface([c3])

        gmsh.model.geo.synchronize()

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Create two materials, one for each side of the fault.
        # We use the `MaterialGroup` data class defined in `gmsh_utils.`
        # The tag argument specifies the integer tag for the physical group.
        # The entities argument specifies the array of surfaces for the material.
        materials = (
            MaterialGroup(tag=1, entities=[self.s_1]),
            MaterialGroup(tag=2, entities=[self.s_2]),
            MaterialGroup(tag=3, entities=[self.s_3]),
            MaterialGroup(tag=4, entities=[self.s_4])
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries and the fault.
        # We use the `VertexGroup` data class defined in `gmsh_utils`.
        # The name and tag specify the name and tag assigned to the physical group.
        # The dimension and entities specify the geometric entities to include in the physical
        # group.
        vertex_groups = (
            VertexGroup(name="boundary_xneg", tag=10, dim=1, entities=[self.c_xneg]),
            VertexGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.c_xpos]),
            VertexGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.c_yneg]),
            VertexGroup(name="boundary_ypos", tag=13, dim=1, entities=[self.c_ypos]),
            VertexGroup(name="fault1", tag=20, dim=1, entities=[self.c_fault1sge1,self.c_fault1sge2,self.c_fault1sge3,self.c_fault1sge4,self.c_fault1sge5,self.c_fault1sge6]),
            VertexGroup(name="fault2", tag=21, dim=1, entities=[self.c_fault2sge1,self.c_fault2sge2,self.c_fault2sge3,self.c_fault2sge4,self.c_fault2sge5,self.c_fault2sge6,self.c_fault2sge7]),
            VertexGroup(name="fault2edge", tag=22, dim=0, entities=[self.p_inter])
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
        field_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_distance, "CurvesList", [self.c_fault1sge1,self.c_fault1sge2,self.c_fault1sge3,
                                                                        self.c_fault1sge4,self.c_fault1sge5,self.c_fault1sge6,
                                                                        self.c_fault2sge1,self.c_fault2sge2,self.c_fault2sge3,
                                                                        self.c_fault2sge4,self.c_fault2sge5,self.c_fault2sge6,self.c_fault2sge7])

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the fault, the distance from
        # the fault (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression` 
        # for creating the string with the mathematical function.
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(field_distance, min_dx=4.0e+3, bias=1.05)
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

