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
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh, group_exclude)

class App(GenerateMesh):
    """
    Application for generating the mesh.
    """
    X_WEST = -600.0e+3
    X_EAST = 600.0e+3
    Y_BOT = -600.0e+3
    Y_MOHO = -40.0e+3

    DX_FAULT = 5.0e+3
    DX_BIAS = 1.07

    # Topography/bathymetry extracted manually from Google Earth
    TOPO_POINTS = (
        (   X_WEST, -2000.0),
        (-439.1e+3, -300.0),
        (-351.2e+3, -800.0),
        (-263.4e+3,     0.0),
        (-175.6e+3,   400.0),
        ( -87.7e+3,     0.0),
        (   0.0e+3,  -400.0),
        (  87.7e+3, -3000.0),
        ( 165.6e+3, -6000.0),
        ( 263.4e+3, -5400.0),
        ( 351.2e+3, -5400.0),
        ( 439.1e+3, -5400.0),
        (   X_EAST, -5700.0),
    )
    TOPO_WEST = 0
    TOPO_TRENCH = 8
    TOPO_EAST = len(TOPO_POINTS)-1

    # Top of slab from Slab 1.0 (Hayes et al., 2012, https://doi.org/10.1029/2011JB008524) 
    SLABTOP_POINTS = (
        (-422.4e+3, -240.00e+3),
        (-331.0e+3, -180.00e+3),
        (-261.6e+3, -140.00e+3),
        (-223.9e+3, -120.00e+3),
        (-182.6e+3, -100.00e+3),
        (-134.3e+3,  -80.00e+3),
        ( -84.5e+3,  -63.00e+3),
        ( -74.6e+3,  -60.00e+3),
        (  -7.9e+3,     Y_MOHO),
        (  71.1e+3,  -20.00e+3),
        ( 160.5e+3,   -7.50e+3),
    )
    SLABTOP_WEST = 0
    SLABTOP_COSEISMIC = 6
    SLABTOP_MOHO = 8
    SLABTOP_EAST = len(SLABTOP_POINTS)-1

    SLABBOT_ADD = (
        (175.6e+3, Y_MOHO),
        (  X_EAST, Y_MOHO),
    )

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
        # Create curve for topography/bathymetry
        pts_topo = []
        for x, y in self.TOPO_POINTS:
            pt = gmsh.model.geo.add_point(x, y, 0.0)
            pts_topo.append(pt)
        c_topo = gmsh.model.geo.add_spline(pts_topo)
        p_topo_west = pts_topo[self.TOPO_WEST]
        p_topo_east = pts_topo[self.TOPO_EAST]
        p_topo_trench = pts_topo[self.TOPO_TRENCH]
        
        # Create b-spline curve for the top of the slab
        pts_slabtop = []
        for x, y in self.SLABTOP_POINTS:
            pt = gmsh.model.geo.add_point(x, y, 0.0)
            pts_slabtop.append(pt)
        c_slabtop = gmsh.model.geo.add_bspline(pts_slabtop + [p_topo_trench])
        self.p_slabtop_west = pts_slabtop[self.SLABTOP_WEST]
        p_slabtop_east = pts_slabtop[self.SLABTOP_EAST]
        self.p_slabtop_coseismic = pts_slabtop[self.SLABTOP_COSEISMIC]
        p_slabtop_moho = pts_slabtop[self.SLABTOP_MOHO]

        # Create b-spline curve for the bottom of the slab.
        # We translate points to the east. A better approach would be to
        # move points normal to the slab to preserve uniform thickness.
        pts_slabbot = []
        offset = 120.0e+3
        for x, y in self.SLABTOP_POINTS[0:8]:
            pt = gmsh.model.geo.add_point(x+offset, y, 0.0)
            pts_slabbot.append(pt)
        for x, y in self.SLABBOT_ADD:
            pt = gmsh.model.geo.add_point(x, y, 0.0)
            pts_slabbot.append(pt)
        self.c_slabbot = gmsh.model.geo.add_bspline(pts_slabbot)
        self.p_slabbot_west = pts_slabbot[0]
        p_slabbot_east = pts_slabbot[-1]

        # Create curve for bottom edge of slab
        c_slab_end = gmsh.model.geo.add_polyline([self.p_slabtop_west, self.p_slabbot_west])


        # Create top of mantle
        p_moho_west = gmsh.model.geo.add_point(self.X_WEST, self.Y_MOHO, 0.0)
        self.c_conmoho = gmsh.model.geo.add_line(p_moho_west, p_slabtop_moho)

        # Create lateral edges and bottom
        p_bot_west = gmsh.model.geo.add_point(self.X_WEST, self.Y_BOT, 0.0)
        p_bot_east = gmsh.model.geo.add_point(self.X_EAST, self.Y_BOT, 0.0)
        c_west = gmsh.model.geo.add_polyline([p_bot_west, p_moho_west, p_topo_west])
        self.c_bot = gmsh.model.geo.add_polyline([p_bot_west, p_bot_east])
        c_east = gmsh.model.geo.add_polyline([p_bot_east, p_slabbot_east, p_topo_east])

        # Split curves to form bounding curves for each material
        # Constructing the entire boundary curves as splines and then breaking them into
        # pieces preserves C1 continuity in the curves.
        curves = gmsh.model.geo.split_curve(c_topo, [p_topo_trench])
        self.c_topo_west = curves[0]
        self.c_topo_east = curves[1]

        curves = gmsh.model.geo.split_curve(c_slabtop, [self.p_slabtop_coseismic, p_slabtop_moho])
        self.c_slabtop_mantle_lower = curves[0]
        self.c_slabtop_mantle_upper = curves[1]
        self.c_slabtop_crust = curves[2]

        curves = gmsh.model.geo.split_curve(c_west, [p_moho_west])
        self.c_west_mantle = curves[0]
        self.c_west_crust = curves[1]

        curves = gmsh.model.geo.split_curve(c_east, [p_slabbot_east])
        self.c_east_mantle = curves[0]
        self.c_east_crust = curves[1]

        # Create surfaces from bounding curves
        loop = gmsh.model.geo.add_curve_loop([
            -self.c_west_crust,
            self.c_conmoho,
            self.c_slabtop_crust,
            -self.c_topo_west,
            ])
        self.s_concrust = gmsh.model.geo.add_plane_surface([loop])

        loop = gmsh.model.geo.add_curve_loop([
            c_slab_end,
            self.c_slabbot,
            self.c_east_crust,
            -self.c_topo_east,
            -self.c_slabtop_crust,
            -self.c_slabtop_mantle_upper,
            -self.c_slabtop_mantle_lower,
            ])
        self.s_oceancrust = gmsh.model.geo.add_plane_surface([loop])

        loop = gmsh.model.geo.add_curve_loop([
            -self.c_west_mantle,
            self.c_bot,
            self.c_east_mantle,
            -self.c_slabbot,
            -c_slab_end,
            self.c_slabtop_mantle_lower,
            self.c_slabtop_mantle_upper,
            -self.c_conmoho,
            ])
        self.s_mantle = gmsh.model.geo.add_plane_surface([loop])


        gmsh.model.geo.synchronize()


    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented.
        """
        # Create materials matching surfaces.
        materials = (
            MaterialGroup(tag=1, entities=[self.s_concrust]),
            MaterialGroup(tag=2, entities=[self.s_oceancrust]),
            MaterialGroup(tag=3, entities=[self.s_mantle]),
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries and the fault.
        vertex_groups = (
            VertexGroup(name="groundsurf", tag=10, dim=1, entities=[self.c_topo_west, self.c_topo_east]),
            VertexGroup(name="bndry_west", tag=11, dim=1, entities=[self.c_west_mantle, self.c_west_crust]),
            VertexGroup(name="bndry_east_crust", tag=12, dim=1, entities=[self.c_east_crust]),
            VertexGroup(name="bndry_east_mantle_tmp", tag=113, dim=1, entities=[self.c_east_mantle]),
            VertexGroup(name="bndry_bot", tag=14, dim=1, entities=[self.c_bot]),
            VertexGroup(name="fault_coseismic", tag=20, dim=1, entities=[self.c_slabtop_mantle_upper, self.c_slabtop_crust]),
            VertexGroup(name="fault_coseismic_edge", tag=30, dim=0, entities=[self.p_slabtop_coseismic]),
            VertexGroup(name="fault_slabtop", tag=21, dim=1, entities=[self.c_slabtop_mantle_lower, self.c_slabtop_mantle_upper, self.c_slabtop_crust]),
            VertexGroup(name="fault_slabtop_edge", tag=31, dim=0, entities=[self.p_slabtop_west]),
            VertexGroup(name="fault_slabbot", tag=22, dim=1, entities=[self.c_slabbot]),
            VertexGroup(name="fault_slabbot_edge", tag=32, dim=0, entities=[self.p_slabbot_west]),
        )
        for group in vertex_groups:
            group.create_physical_group()
        group_exclude("bndry_east_mantle_tmp", "fault_slabbot", new_name="bndry_east_mantle", new_tag=13)

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
        gmsh.model.mesh.field.setNumbers(fault_distance, "CurvesList", [self.c_slabtop_mantle_upper, self.c_slabtop_crust])

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the fault, the distance from
        # the fault (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression` 
        # for creating the string with the mathematical function.
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(fault_distance, min_dx=self.DX_FAULT, bias=self.DX_BIAS)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        # Finally, we use the field `field_size` for the cell size of the mesh.
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
