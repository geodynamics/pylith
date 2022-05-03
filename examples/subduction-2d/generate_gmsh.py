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
    Y_BOT = -340.0e+3
    Y_MOHO = -40.0e+3

    # Topography/bathymetry extracted manually from Google Earth
    TOPOBATHY_POINTS = (
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
    TOPOBATHY_WEST = 0
    TOPOBATHY_TRENCH = 8
    TOPOBATHY_EAST = len(TOPOBATHY_POINTS)-1

    # Top of slab from Slab 1.0 (Hayes et al., 2012, https://doi.org/10.1029/2011JB008524) 
    SLABTOP_POINTS = (
        (   X_WEST,      Y_BOT),
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
    SLABTOP_COSEISMIC = 7
    SLABTOP_MOHO = 9
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

    def create_geometry(self):
        """Create geometry.
        """
        # Create curve for topography/bathymetry
        pts_topobathy = []
        for x, y in self.TOPOBATHY_POINTS:
            pt = gmsh.model.geo.add_point(x, y, 0.0)
            pts_topobathy.append(pt)
        l_topobathy = gmsh.model.geo.add_spline(pts_topobathy)
        p_topobathy_west = pts_topobathy[self.TOPOBATHY_WEST]
        p_topobathy_east = pts_topobathy[self.TOPOBATHY_EAST]
        p_topobathy_trench = pts_topobathy[self.TOPOBATHY_TRENCH]
        
        # Create b-spline curve for the top of the slab
        pts_slabtop = []
        for x, y in self.SLABTOP_POINTS:
            pt = gmsh.model.geo.add_point(x, y, 0.0)
            pts_slabtop.append(pt)
        l_slabtop = gmsh.model.geo.add_bspline(pts_slabtop + [p_topobathy_trench])
        p_slabtop_west = pts_slabtop[self.SLABTOP_WEST]
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
        self.l_slabbot = gmsh.model.geo.add_bspline(pts_slabbot)
        p_slabbot_west = pts_slabbot[0]
        p_slabbot_east = pts_slabbot[-1]

        # Create top of mantle
        p_moho_west = gmsh.model.geo.add_point(-600.0e+3, -40.0e+3, 0.0)
        self.l_conmoho = gmsh.model.geo.add_line(p_moho_west, p_slabtop_moho)

        # Create lateral edges and bottom
        p_bot_east = gmsh.model.geo.add_point(600.0e+3, -340.0e+3, 0.0)
        l_west = gmsh.model.geo.add_polyline([p_slabtop_west, p_moho_west, p_topobathy_west])
        l_bot = gmsh.model.geo.add_polyline([p_slabtop_west, p_slabbot_west, p_bot_east])
        l_east = gmsh.model.geo.add_polyline([p_bot_east, p_slabbot_east, p_topobathy_east])

        # Split curves to form bounding curves for each material
        # Constructing the entire boundary curves as splines and then breaking them into
        # pieces preserves C1 continuity in the curves.
        curves = gmsh.model.geo.split_curve(l_topobathy, [p_topobathy_trench])
        self.l_topography_west = curves[0]
        self.l_topography_east = curves[1]

        curves = gmsh.model.geo.split_curve(l_slabtop, [self.p_slabtop_coseismic, p_slabtop_moho])
        self.l_slabtop_mantle_lower = curves[0]
        self.l_slabtop_mantle_upper = curves[1]
        self.l_slabtop_crust = curves[2]

        curves = gmsh.model.geo.split_curve(l_west, [p_moho_west])
        self.l_west_mantle = curves[0]
        self.l_west_crust = curves[1]

        curves = gmsh.model.geo.split_curve(l_bot, [p_slabbot_west])
        self.l_bot_slab = curves[0]
        self.l_bot_mantle = curves[1]

        curves = gmsh.model.geo.split_curve(l_east, [p_slabbot_east])
        self.l_east_mantle = curves[0]
        self.l_east_crust = curves[1]

        # Create surfaces from bounding curves
        loop = gmsh.model.geo.add_curve_loop([
            -self.l_west_crust,
            self.l_conmoho,
            self.l_slabtop_crust,
            -self.l_topography_west,
            ])
        self.s_concrust = gmsh.model.geo.add_plane_surface([loop])

        loop = gmsh.model.geo.add_curve_loop([
            -self.l_west_mantle,
            self.l_slabtop_mantle_lower,
            self.l_slabtop_mantle_upper,
            -self.l_conmoho,
            ])
        self.s_conmantle = gmsh.model.geo.add_plane_surface([loop])

        loop = gmsh.model.geo.add_curve_loop([
            self.l_bot_slab,
            self.l_slabbot,
            self.l_east_crust,
            -self.l_topography_east,
            -self.l_slabtop_crust,
            -self.l_slabtop_mantle_upper,
            -self.l_slabtop_mantle_lower,
            ])
        self.s_oceancrust = gmsh.model.geo.add_plane_surface([loop])

        loop = gmsh.model.geo.add_curve_loop([
            self.l_bot_mantle,
            self.l_east_mantle,
            -self.l_slabbot,
            ])
        self.s_oceanmantle = gmsh.model.geo.add_plane_surface([loop])

        gmsh.model.geo.synchronize()


    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented.
        """
        # Create materials matching surfaces.
        materials = (
            MaterialGroup(tag=1, entities=[self.s_concrust]),
            MaterialGroup(tag=2, entities=[self.s_conmantle]),
            MaterialGroup(tag=3, entities=[self.s_oceancrust]),
            MaterialGroup(tag=4, entities=[self.s_oceanmantle]),
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries and the fault.
        vertex_groups = (
            VertexGroup(name="groundsurf", tag=10, dim=1, entities=[self.l_topography_west, self.l_topography_east]),
            VertexGroup(name="bndry_west_incslab", tag=12, dim=1, entities=[self.l_west_mantle, self.l_west_crust]),
            VertexGroup(name="bndry_east_crust", tag=13, dim=1, entities=[self.l_east_crust]),
            VertexGroup(name="bndry_east_mantle", tag=14, dim=1, entities=[self.l_east_mantle]),
            VertexGroup(name="bndry_bot_slab", tag=15, dim=1, entities=[self.l_bot_slab]),
            VertexGroup(name="bndry_bot_mantle", tag=16, dim=1, entities=[self.l_bot_mantle]),
            VertexGroup(name="fault_coseismic", tag=20, dim=1, entities=[self.l_slabtop_mantle_upper, self.l_slabtop_crust]),
            VertexGroup(name="fault_coseismic_edge", tag=21, dim=0, entities=[self.p_slabtop_coseismic]),
            VertexGroup(name="fault_slabtop", tag=22, dim=1, entities=[self.l_slabtop_mantle_lower, self.l_slabtop_mantle_upper, self.l_slabtop_crust]),
            VertexGroup(name="fault_slabbot", tag=23, dim=1, entities=[self.l_slabbot]),
        )
        for group in vertex_groups:
            group.create_physical_group()
        group_exclude("bndry_west_incslab", "fault_slabtop", new_name="bndry_west", new_tag=11)

    def generate_mesh(self, cell):
        """Generate the mesh.
        """
        # Set discretization size with geometric progression from distance to the fault

        # First, we setup a field `field_distance` with the distance from the fault.
        fault_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(fault_distance, "CurvesList", [self.l_slabtop_mantle_upper, self.l_slabtop_crust])

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the fault, the distance from
        # the fault (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression` 
        # for creating the string with the mathematical function.
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(fault_distance, min_dx=5.0e+3, bias=1.07)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        ## Finally, we use the field `field_size` for the cell size of the mesh.
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        gmsh.option.setNumber("Mesh.Algorithm", 8)
        gmsh.model.mesh.generate(2)
        if cell == "quad":
            gmsh.model.mesh.recombine()
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()


# End of file
