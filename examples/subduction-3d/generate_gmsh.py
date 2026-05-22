import numpy as np
import gzip
# Import Gmsh Python interface
import gmsh

# Import the gmsh_utils Python module supplied with PyLith.
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)
from pyproj import Transformer, CRS
import netCDF4
import math


class App(GenerateMesh):
    """
    Application used to generate the mesh using Gmsh.

    App uses `GenerateMesh` from `gmsh_utils` for common functionality that we avoid
    duplicating in each of our examples.
    """
    CONTOURS_FILENAME = "cas_contours_dep.in.txt.gz"
    SLAB_THICKNESS = 50.0 * 1000
    FILENAME_LOCALDEM = "topography.nc"
    BOX_SIDE_LENGTH = 800*1000
    BOX_DEPTH_LENGTH = 400*1000

    UP_DIP_ELEV = 6.0*1000
    UP_DIP_DIST = 600.0*1000
    UP_DIP_ANGLE = math.radians(10.0)
    FAULT_STRIKE = math.radians(0.0)

    DX_MIN = 1.0e+4
    DX_BIAS = 1.1

    CRUST_DEPTH = 40.0 * 1000
    PATCH_LENGTH = 200.0 * 1000

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

    def _read_and_generate_splines(self):
        """Read geographic coordinates from "cas_contours_dep.in.txt.gz" and return a dictionary where the keys are
        the depth of the contour and the values are the coordinates (latitude, longitude, elevation) as a numpy array.
        """
        with gzip.open(self.CONTOURS_FILENAME, "rb") as file:
            lines = file.readlines()
        contours = {}
        points = []
        key = None
        all_points = []
        for line in lines:
            if line.decode().strip() == "END":
                contours[key] = np.array(points, dtype=np.float64)
                points = []
                continue
            if len(line.split()) == 1:
                key = int(line)
                continue
            pt = list(map(float, line.strip().split()))  # lon/lat/elev
            this_point = [pt[1], pt[0], pt[2]*1000]
            points.append(this_point)  # lat/lon/elev
            all_points.append(this_point)
        return contours, np.array(all_points)

    def _generate_extended_contours(self, contours):
        """Add contours up-dip from original contours.

        We increase the horizontal distance between the contours at a
        geometric rate. The first contour is at a distance of
        dist_horiz, followed by 2*dist_horiz, 4*dist_horiz, etc.

        The horizontal distance of contour n from the original one is
        (2**(n+1)-1)) * dist_horiz, n=0,1,2,...
        """
        key = min(contours.keys())
        contour_top = contours[key]
        z_top = contour_top[0][2]
        dist_horiz = (self.UP_DIP_ELEV - z_top) / math.tan(self.UP_DIP_ANGLE)
        dx = -dist_horiz * math.cos(self.FAULT_STRIKE)
        dy = dist_horiz * math.sin(self.FAULT_STRIKE)
        contours_up_dip = {}
        ncontours = int(math.ceil(math.log((self.UP_DIP_DIST / dist_horiz) + 1) / math.log(2.0)))
        target_range = np.linspace(z_top, self.UP_DIP_ELEV, ncontours+1)
        for i in range(ncontours):
            contour = np.array(contour_top)
            contour[:, 0] += (2 ** i) * dx
            contour[:, 1] += (2 ** i) * dy
            contour[:, 2] = target_range[i+1]
            contours_up_dip[-target_range[i+1]/1000] = contour
        return contours_up_dip

    @staticmethod
    def _generate_gmsh_contour(contour: np.ndarray):
        """Helper function for creating a spline and wire from a np array of shape (x,3)
        """
        points = []
        for point in contour:
            gmsh_point = gmsh.model.occ.add_point(point[0], point[1], point[2])
            points.append(gmsh_point)
        spline = gmsh.model.occ.add_spline(points)
        wire = gmsh.model.occ.add_wire([spline])
        return spline, wire

    @staticmethod
    def add_plane_surface_at_point(x, y, z, x_extent, y_extent):
        """Helper function for creating a plane surface
        """
        p1 = gmsh.model.occ.add_point(x, y, z)
        p2 = gmsh.model.occ.add_point(x + x_extent, y, z)
        p3 = gmsh.model.occ.add_point(x + x_extent, y + y_extent, z)
        p4 = gmsh.model.occ.add_point(x, y + y_extent, z)

        l1 = gmsh.model.occ.add_line(p1, p2)
        l2 = gmsh.model.occ.add_line(p2, p3)
        l3 = gmsh.model.occ.add_line(p3, p4)
        l4 = gmsh.model.occ.add_line(p4, p1)

        loop = gmsh.model.occ.add_curve_loop([l1, l2, l3, l4])
        surface = gmsh.model.occ.add_plane_surface([loop])
        return surface

    @staticmethod
    def _calculate_bounding_box(lats, longs, radius=0):
        """Helper function for calculating the bounding box, takes in an array of shape (x,3)
        """
        # Calculate the min and max for latitude and longitude
        min_lat = np.min(lats) - radius
        max_lat = np.max(lats) + radius
        min_lon = np.min(longs) - radius
        max_lon = np.max(longs) + radius

        return min_lon, max_lat, max_lon, min_lat

    def create_geometry(self):
        """Create geometry.
        We use a local Transverse Mercator projection centered at latitude 45.5231 and longitude -122.6765.
        We use contour traces in geographic coordinates from "cas_contours_dep.in.txt.gz"
        """
        crs_wgs84_3d = CRS.from_epsg(4979)
        crs_projected = CRS.from_proj4(
            "+proj=tmerc +datum=WGS84 +lat_0=45.5231 +lon_0=-122.6765 +k=0.9996 +units=m +type=crs"
        )
        transformer = Transformer.from_crs(crs_wgs84_3d, crs_projected, always_xy=True)

        dem_nc = netCDF4.Dataset(self.FILENAME_LOCALDEM)
        latitude = dem_nc.variables["lat"][:]
        longitude = dem_nc.variables["lon"][:]
        elevation = dem_nc.variables["Band1"][:].astype(float)

        # Generate topography surface
        topo_mask = elevation.recordmask
        wires = []
        for i in range(elevation.shape[0]):
            points = []
            for j in range(elevation.shape[1]):
                if isinstance(topo_mask, np.ndarray):
                    if topo_mask[i, j]:
                        continue  # Skip masked values
                point_latitude = latitude[i]
                point_longitude = longitude[j]
                point_elevation = elevation[i, j]
                position = transformer.transform(point_longitude, point_latitude, point_elevation)
                points.append(position)
            points = np.array(points)
            _, wire = self._generate_gmsh_contour(points)
            wires.append(wire)
        topography_surface = gmsh.model.occ.add_thru_sections(wires, makeSolid=False, makeRuled=False, maxDegree=2)

        # Read slab contours
        contours, all_points = self._read_and_generate_splines()
        bounding_box = self._calculate_bounding_box(lats=all_points[:, 0], longs=all_points[:, 1])
        print("Bounding box:", bounding_box)
        contours[5] = np.flip(contours[5], axis=0)  # The orientation of the first contour in this file is reversed

        # Project to mercator
        projected_contours = {}
        for depth, contour in contours.items():
            position = transformer.transform(contour[:, 1], contour[:, 0], contour[:, 2])
            projected_contours[depth] = np.array(position).T

        # Extend the slab to the topography
        projected_contours_up_dip = self._generate_extended_contours(projected_contours)
        projected_contours = projected_contours_up_dip | projected_contours  # merge the 2 dicts

        # Generate the top surface of the slab
        all_points_on_grid = []
        first_side = []
        third_side = []
        sorted_project_contours = [v for k, v in sorted(projected_contours.items())]
        for projected_contour in sorted_project_contours:
            gmsh_points = []
            for point in projected_contour:
                gmsh_point = gmsh.model.occ.add_point(point[0], point[1], point[2])
                gmsh_points.append(gmsh_point)
            first_side.append(gmsh_points[0])
            third_side.append(gmsh_points[-1])
            all_points_on_grid.append(gmsh_points)

        first_side_spline = gmsh.model.occ.add_spline(first_side)
        third_side_spline = gmsh.model.occ.add_spline(third_side)
        second_side_spline = gmsh.model.occ.add_spline(all_points_on_grid[0])
        forth_side_spline = gmsh.model.occ.add_spline(all_points_on_grid[-1])

        wire = gmsh.model.occ.add_wire([first_side_spline, second_side_spline, third_side_spline, forth_side_spline])
        slab_top_surface = gmsh.model.occ.addSurfaceFilling(wire, pointTags=np.concatenate(all_points_on_grid).tolist())

        gmsh.model.occ.synchronize()

        # Find the normals of the slab surface
        space = np.linspace(0, 1, 10)
        xv, yv = np.meshgrid(space, space)
        parametricCoords = np.column_stack((xv.flatten(), yv.flatten())).ravel()
        normals = gmsh.model.getNormal(slab_top_surface, parametricCoords)
        normal = np.average(normals.reshape(-1, 3), axis=0)

        # Extrude the slab surface to get it's volume
        dx = self.SLAB_THICKNESS * normal[0]
        dy = self.SLAB_THICKNESS * normal[1]
        dz = self.SLAB_THICKNESS * normal[2]
        slab_volume = gmsh.model.occ.extrude([(2, slab_top_surface)], dx=dx, dy=dy, dz=dz, recombine=True)

        # Generate the splay fault
        splay_bottom = np.copy(projected_contours[15])
        splay_top = np.copy(projected_contours[15])

        splay_bottom[:, 2] -= 8.0e+3
        _, splay_bottom_wire = self._generate_gmsh_contour(splay_bottom)

        splay_top[:, 2] = 3.0e+3
        splay_top[:, 0] -= 24.0e+3
        _, splay_top_wire = self._generate_gmsh_contour(splay_top)
        splay_surface = gmsh.model.occ.add_thru_sections(
            [splay_bottom_wire, splay_top_wire], makeSolid=False, makeRuled=False, maxDegree=2)

        # Generate the bounding box
        bounding_box = gmsh.model.occ.add_box(
            -self.BOX_SIDE_LENGTH/2,
            -self.BOX_SIDE_LENGTH/2,
            30*1000-self.BOX_DEPTH_LENGTH,
            self.BOX_SIDE_LENGTH,
            self.BOX_SIDE_LENGTH,
            self.BOX_DEPTH_LENGTH
        )

        # Create the topography on top of the bounding box
        gmsh.model.occ.fragment(topography_surface, [[3, bounding_box]])
        gmsh.model.occ.remove([(3, 3), (2, 20)], recursive=True)

        # Fragment the bounding box with the slab to get the main volume and the slab volume inside the bounding box
        gmsh.model.occ.fragment([(3, 2)], [(3, 1)])
        gmsh.model.occ.remove([(3, 3)], recursive=True)

        # Fragment the main volume with the splay surface to get the wedge volume
        gmsh.model.occ.fragment([(3, 1)], splay_surface)

        # Create the crust surface
        crust_surface = self.add_plane_surface_at_point(
            -self.BOX_SIDE_LENGTH / 2 - 50 * 1000,
            -self.BOX_SIDE_LENGTH / 2 - 50 * 1000,
            -1 * 1000 - self.CRUST_DEPTH,  # Crust side of the bounding box about is approximately -1km below 0
            self.BOX_SIDE_LENGTH + 100 * 1000,
            self.BOX_SIDE_LENGTH + 100 * 1000,
        )

        # Split the main volume to get the crust volume
        gmsh.model.occ.fragment([(3, 3)], [(2, crust_surface)])
        gmsh.model.occ.remove_all_duplicates()
        gmsh.model.occ.synchronize()

        # Generate patch
        patch_box = gmsh.model.occ.add_box(
            -self.BOX_SIDE_LENGTH / 2,
            -self.PATCH_LENGTH / 2,
            -400.0 * 1000,
            1000.0 * 1000,
            self.PATCH_LENGTH,
            600.0 * 1000
        )

        gmsh.model.occ.remove([(3, patch_box)])
        gmsh.model.occ.remove([(2, 54), (2, 55), (2, 58), (2, 59)], recursive=True)
        gmsh.model.occ.fragment([(3, 2)], [(2, 56), (2, 57)])
        gmsh.model.occ.remove_all_duplicates()

        gmsh.model.occ.synchronize()
        # gmsh.fltk.run()

        # Get tags of the different sections using GUI
        self.mantle_volume = [5]
        self.crust_volume = [6]
        self.wedge_volume = [4]
        self.slab_volume = [7, 8, 9]

        self.surface_west_xpos = [86, 89]
        self.surface_south_yneg = [88, 77, 83, 51]
        self.surface_south_no_slab = [88, 77, 83]
        self.surface_top_zpos = [52, 60, 68] + [78] + [90]
        self.surface_north_ypos = [81, 84, 91, 70]
        self.surface_north_no_slab = [81, 84, 91]
        self.surface_bottom_zneg = [85]
        self.surface_east_xneg = [80, 82] + [50, 59, 67]
        self.surface_east_no_slab = [80, 82]

        self.slab_top = [52, 56, 57, 58] + [60, 66, 65, 64] + [68, 69, 74, 73]
        self.slab_top_edge = [177, 191, 200]
        self.slab_bottom = [54, 62, 71]
        self.slab_bottom_edge = [176, 190, 199]
        self.slab_top_patch = [66, 65]
        self.slab_top_patch_edge = [175, 174, 192, 188, 189]

        self.surface_splay = [56, 66, 69]
        self.edge_splay = [178, 193, 198]

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.
        """
        materials = (
            MaterialGroup(tag=1, entities=self.slab_volume),
            MaterialGroup(tag=2, entities=self.wedge_volume),
            MaterialGroup(tag=3, entities=self.mantle_volume),
            MaterialGroup(tag=4, entities=self.crust_volume),
        )

        for material in materials:
            material.create_physical_group()

        face_groups = (
            BoundaryGroup(name="boundary_yneg", tag=10, dim=2, entities=self.surface_south_yneg),
            BoundaryGroup(name="boundary_xneg", tag=11, dim=2, entities=self.surface_east_xneg),
            BoundaryGroup(name="boundary_ypos", tag=12, dim=2, entities=self.surface_north_ypos),
            BoundaryGroup(name="boundary_xpos", tag=13, dim=2, entities=self.surface_west_xpos),
            BoundaryGroup(name="boundary_zneg", tag=14, dim=2, entities=self.surface_bottom_zneg),
            BoundaryGroup(name="boundary_zpos", tag=15, dim=2, entities=self.surface_top_zpos),
            BoundaryGroup(name="boundary_xneg_noslab", tag=16, dim=2, entities=self.surface_east_no_slab),
            BoundaryGroup(name="boundary_yneg_noslab", tag=17, dim=2, entities=self.surface_south_no_slab),
            BoundaryGroup(name="boundary_ypos_noslab", tag=18, dim=2, entities=self.surface_north_no_slab),

            BoundaryGroup(name="fault_slabtop", tag=20, dim=2, entities=self.slab_top),
            BoundaryGroup(name="fault_slabbot", tag=21, dim=2, entities=self.slab_bottom),
            BoundaryGroup(name="fault_slabtop_edge", tag=22, dim=1, entities=self.slab_top_edge),
            BoundaryGroup(name="fault_slabbot_edge", tag=23, dim=1, entities=self.slab_bottom_edge),
            BoundaryGroup(name="fault_slabtop_patch", tag=24, dim=2, entities=self.slab_top_patch),
            BoundaryGroup(name="fault_slabtop_patch_edge", tag=25, dim=1, entities=self.slab_top_patch_edge),

            BoundaryGroup(name="fault_splay", tag=30, dim=2, entities=self.surface_splay),
            BoundaryGroup(name="fault_splay_edge", tag=31, dim=1, entities=self.edge_splay),
        )
        for group in face_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh.
        """
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # We setup a field `field_distance` with the distance from the top of the slab.
        field_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_distance, "SurfacesList", self.slab_top)

        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(field_distance, min_dx=self.DX_MIN, bias=self.DX_BIAS)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Laplace2D")


# If script is called from the command line, run the application.
if __name__ == "__main__":
    App().main()

# End of file
