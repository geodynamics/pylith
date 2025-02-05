#!/usr/bin/env nemesis
"""Generate a tri or quad mesh of a subduction zone vertical profile using Gmsh, making
use of the built-in geometry engine.

Points have been projected from longitude/latitude into a local
transverse Mercator projection. PyLith uses the Proj.4 library
for geographic projections. The proj parameters are:

+proj=tmerc +datum=WGS84 +lon_0=142.0 +lat_0=38.0 +k=0.9996

so that the local origin is at a longitude of 142.0 degrees (WGS84)
and a latitude of 38.0 degrees (WGS84).

Run `generate_gmsh.py --help` to see the command line options.
"""
import gmsh
import numpy as np
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Application for generating the mesh.
    """
    # The elastic thickness of the Aleutian Plate is 29.7 +/- 5km (Watts & Zhong, 2002).
    X_WEST = 0.0e+3
    X_EAST = 150.0e+3
    Y_BOT = -30.0e+3
    Y_TOP = 0.0e+3

    fault_dip = np.deg2rad(60)
    max_fault_depth = Y_BOT / 2

    OUTER_RISE_FAULT_SURFACES_X = np.array([100.0e+3, 110.0e+3, 120.0e+3])
    OUTER_RISE_FAULT_SURFACES_Y = np.zeros(len(OUTER_RISE_FAULT_SURFACES_X))

    OUTER_RISE_FAULT_BURIED_EDGES_Y = np.ones(len(OUTER_RISE_FAULT_SURFACES_X)) * max_fault_depth
    OUTER_RISE_FAULT_BURIED_EDGES_X = OUTER_RISE_FAULT_SURFACES_X - OUTER_RISE_FAULT_BURIED_EDGES_Y / np.tan(fault_dip)

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
        """
        p_LL = gmsh.model.geo.add_point(self.X_WEST, self.Y_BOT, 0.0)
        p_LR = gmsh.model.geo.add_point(self.X_EAST, self.Y_BOT, 0.0)
        p_UR = gmsh.model.geo.add_point(self.X_EAST, self.Y_TOP, 0.0)
        p_UL = gmsh.model.geo.add_point(self.X_WEST, self.Y_TOP, 0.0)

        # Create Domain Boundary Curves
        self.c_left = gmsh.model.geo.add_line(p_UL, p_LL)
        self.c_bottom = gmsh.model.geo.add_line(p_LL, p_LR)
        self.c_right = gmsh.model.geo.add_line(p_LR, p_UR)

        self.c_outer_rise_faults = np.zeros(len(self.OUTER_RISE_FAULT_SURFACES_X), dtype=int)
        self.outer_rise_fault_surface_points = np.zeros(len(self.OUTER_RISE_FAULT_SURFACES_X), dtype=int)
        self.outer_rise_fault_buried_edges = np.zeros(len(self.OUTER_RISE_FAULT_BURIED_EDGES_X), dtype=int)

        for i in range(len(self.OUTER_RISE_FAULT_SURFACES_X)):
            self.outer_rise_fault_buried_edges[i] = gmsh.model.geo.add_point(self.OUTER_RISE_FAULT_BURIED_EDGES_X[i],
                                                                             self.OUTER_RISE_FAULT_BURIED_EDGES_Y[i],
                                                                             0.0)

            self.outer_rise_fault_surface_points[i] = gmsh.model.geo.add_point(self.OUTER_RISE_FAULT_SURFACES_X[i],
                                                                               self.OUTER_RISE_FAULT_SURFACES_Y[i],
                                                                               0.0)

            self.c_outer_rise_faults[i] = gmsh.model.geo.add_polyline([self.outer_rise_fault_surface_points[i], self.outer_rise_fault_buried_edges[i]])
        
        self.subducting_top_bound_points = p_UR
        self.subducting_top_bound_points = np.insert(self.subducting_top_bound_points, 0, self.outer_rise_fault_surface_points)
        self.subducting_top_bound_points = np.insert(self.subducting_top_bound_points, 0, p_UL)

        self.c_top = gmsh.model.geo.add_polyline(self.subducting_top_bound_points)

        self.all_curves = np.zeros(len(self.outer_rise_fault_surface_points) + 1, dtype=int)
        curves = gmsh.model.geo.split_curve(self.c_top, [self.outer_rise_fault_surface_points[0]])
        self.all_curves[0] = curves[0]
        self.all_curves[1] = curves[1]
        for i in range(1, len(self.outer_rise_fault_surface_points)):
            curves = gmsh.model.geo.split_curve(self.all_curves[i], [self.outer_rise_fault_surface_points[i]])
            self.all_curves[i] = curves[0]
            self.all_curves[i + 1] = curves[1]

        # Create surfaces from bounding curves
        loop_array = np.array([self.c_left,
                               self.c_bottom,
                               self.c_right])

        for j in range(len(self.outer_rise_fault_surface_points)):
            loop_array = np.append(loop_array, -self.all_curves[j])
            loop_array = np.append(loop_array, self.c_outer_rise_faults[j])
            loop_array = np.append(loop_array, -self.c_outer_rise_faults[j])
        loop_array = np.append(loop_array, -self.all_curves[-1])

        loop = gmsh.model.geo.add_curve_loop(loop_array)
        self.s_slab = gmsh.model.geo.add_plane_surface([loop])
        gmsh.model.geo.synchronize()


    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented.
        """
        # Create materials matching surfaces.
        materials = (
            MaterialGroup(tag=1, entities=[self.s_slab]),
        )
        for material in materials:
            material.create_physical_group()

        top_boundary_entities = self.all_curves

        # Create physical groups for the boundaries and the fault.
        face_groups = (
            BoundaryGroup(name="boundary_ypos", tag=10, dim=1, entities=top_boundary_entities),
            BoundaryGroup(name="boundary_xneg", tag=11, dim=1, entities=[self.c_left]),
            BoundaryGroup(name="boundary_xpos", tag=12, dim=1, entities=[self.c_right]),
            BoundaryGroup(name="boundary_yneg", tag=13, dim=1, entities=[self.c_bottom]),

        )
        for i in range(len(self.outer_rise_fault_surface_points)):
            face_groups = np.append(face_groups, BoundaryGroup(name="fault_" + str(i), \
                                                                 tag=201 + i, dim=1, entities=[self.c_outer_rise_faults[i]]))

            face_groups = np.append(face_groups, BoundaryGroup(name="edge_fault_" + str(i), tag=301 + i, dim=0, \
                                                                 entities=[int(self.outer_rise_fault_buried_edges[i])]))
        for group in face_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh.
        """
        # Set discretization size with geometric progression from distance to the fault.
        # We turn off the default sizing methods.
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # First, we setup a field `field_distance` with the distance from the fault.
        fault_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumber(fault_distance, "Sampling", 200)

        mesh_refinement_array = np.empty(0)
        for fault in self.c_outer_rise_faults:
            mesh_refinement_array = np.append(mesh_refinement_array, fault)
        for curve in self.all_curves:
            mesh_refinement_array = np.append(mesh_refinement_array, curve)
        gmsh.model.mesh.field.setNumbers(fault_distance, "CurvesList", mesh_refinement_array)

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the fault, the distance from
        # the fault (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression`
        # for creating the string with the mathematical function.

        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(fault_distance, min_dx=0.5e+3, bias=1.05)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        ## Finally, we use the field `field_size` for the cell size of the mesh.
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        if cell == "quad":
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.recombine()
        else:
            gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")

if __name__ == "__main__":
    App().main()

# End of file
