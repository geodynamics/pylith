#!/usr/bin/env nemesis

import math
import gzip

import numpy
import netCDF4
import gmsh

from pylith.meshio.gmsh_utils import BoundaryGroup, MaterialGroup, GenerateMesh

# Local Python module with coordinate transformations
import coordsys


class App(GenerateMesh):
    """
    Application used to generate the mesh using Gmsh.

    App uses `GenerateMesh` from `gmsh_utils` for common functionality to avoid
    duplicate code in our examples.
    """

    KM = 1000.0  # scale to convert km to m

    SLAB_CONTOURS_FILENAME = "cas_contours_dep.in.txt.gz"
    FILENAME_LOCALDEM = "etopo2020_bedrock_local.nc"
    SLAB_THICKNESS = 50.0 * KM
    DOMAIN_X = DOMAIN_Y = 800 * KM
    DOMAIN_Z = 400 * KM
    TOP_BOX = 10.0 * KM

    # Elevation and dip angle of extra up-dip contours
    UP_DIP_ELEV = 12.0 * KM
    UP_DIP_ANGLE = math.radians(20.0)

    # Discretization size
    DX_MIN = 20.0 * KM
    DX_BIAS = 1.08

    # Depth of bottom of continental crust
    CRUST_DEPTH = 30.0 * KM

    # Size of patch on slab interface where we prescribe coseismic slip
    PATCH_LENGTH = 200.0 * KM
    PATCH_WIDTH = 275.0 * KM

    def __init__(self):
        """Constructor."""
        super().__init__()
        self.debug = True

        # Set the cell choices available through command line options.
        self.cell_choices = {
            "default": "tet",
            "choices": ["tet"],
        }

    def create_geometry(self):
        """Create geometry.
        We use a local Transverse Mercator projection centered at latitude 45.5231 and longitude -122.6765.
        We use contour traces in geographic coordinates from "cas_contours_dep.in.txt.gz"
        """
        topography_surface = self._create_topography_surface()
        slab_top_surface, contours = self._create_slab_top_surface()
        gmsh.model.occ.synchronize()  # needed for surface normal

        slab_normal = self._calculate_surface_avg_normal(slab_top_surface)
        slab_volume = self._create_slab_volume(slab_top_surface, slab_normal)
        splay_fault_surface = self._create_splay_fault_surface(contours)
        domain_box = self._create_domain_box()
        continental_moho = self._create_continental_moho_surface()

        gmsh.model.occ.synchronize()

        # Create the topography on top of the bounding box
        tags, _ = gmsh.model.occ.fragment(topography_surface, [(3, domain_box)])
        domain_tag = tags[0]
        if self.debug:
            print(
                f"Fragment domain with topography, tags={tags}. Setting domain={domain_tag}, deleting={tags[1:]}"
            )
        gmsh.model.occ.remove(tags[1:], recursive=True)

        # Fragment the bounding box with the slab to get the main volume and the slab volume inside the bounding box
        tags, _ = gmsh.model.occ.fragment([domain_tag], [slab_volume])
        if self.debug:
            print("Fragment domain with slab.")
        domain_tag = tags[0]
        slab_tag = tags[1]
        if self.debug:
            print(
                f"Fragment domain with slab, tags={tags}. Setting domain={domain_tag}, slab={slab_tag}, deleting={tags[2:]}"
            )
        gmsh.model.occ.remove(tags[2:], recursive=True)

        # Fragment the main volume with the splay surface to get the wedge volume
        tags, _ = gmsh.model.occ.fragment([domain_tag], splay_fault_surface)
        domain_tag = tags[0]
        wedge_tag = tags[1]
        if self.debug:
            print(
                f"Fragment domain with splay fault, tags={tags}. Setting domain={domain_tag}, slab={slab_tag} wedge={wedge_tag}, deleting={tags[2:]}"
            )
        gmsh.model.occ.remove(tags[2:], recursive=True)

        # Split the main volume to get the crust volume
        tags, _ = gmsh.model.occ.fragment([domain_tag], [(2, continental_moho)])
        domain_tag = tags[0]
        crust_tag = tags[1]
        if self.debug:
            print(
                f"Fragment domain with moho, tags={tags}. Setting domain={domain_tag}, crust={crust_tag}, slab={slab_tag} wedge={wedge_tag}, deleting={tags[2:]}"
            )
        gmsh.model.occ.remove(tags[2:], recursive=True)

        patch_box = gmsh.model.occ.add_box(
            -self.DOMAIN_X / 2,
            -self.PATCH_LENGTH / 2,
            -100.0 * self.KM,
            self.PATCH_WIDTH,
            self.PATCH_LENGTH,
            200.0 * self.KM,
        )

        tags, _ = gmsh.model.occ.fragment([slab_tag], [(3, patch_box)])
        slab_tags = tags[0:2]
        if self.debug:
            print(
                f"Fragment slab with patch, tags={tags}. Setting domain={domain_tag}, crust={crust_tag}, slab={slab_tags} wedge={wedge_tag}, deleting={tags[2:]}"
            )
        gmsh.model.occ.remove(tags[2:], recursive=True)

        # Imprint all volumes on other volumes and remove duplicate entities
        gmsh.model.occ.remove_all_duplicates()
        gmsh.model.occ.synchronize()

        # Get tags of the different sections using GUI

        # Volumes
        self.mantle_volume = [domain_tag[1]]
        self.crust_volume = [crust_tag[1]]
        self.wedge_volume = [wedge_tag[1]]
        self.slab_volume = [v for d, v in slab_tags]

        # Boundary surfaces
        self.surface_east = [78, 81]
        self.surface_west = [74, 86, 90, 94]
        self.surface_west_no_slab = [74]

        self.surface_north = [73, 76, 85, 92]
        self.surface_north_no_slab = [73, 76, 85]
        self.surface_north_no_slab_splay = [76, 85]
        self.surface_south = [67, 75, 80, 87]
        self.surface_south_no_slab = [67, 75, 80]
        self.surface_south_no_slab_splay = [75, 80]

        self.surface_top = [68, 82, 88, 93, 95]
        self.surface_bottom = [77]

        # Faults
        self.slab_top = [69, 71, 79, 83] + [70, 84]
        # self.slab_top_edge = [118] # mac
        self.slab_top_edge = [96]  # linux

        self.slab_bottom = [57, 66]
        # self.slab_bottom_edge = [109] # mac
        self.slab_bottom_edge = [87]  # linux

        self.slab_top_patch = [70, 84]
        # self.slab_top_patch_edge = [200, 223, 182, 224, 203] # mac
        self.slab_top_patch_edge = [160, 178, 181, 201, 202]  # linux

        self.surface_splay = [72]
        # self.splay_edge = [201, 204, 205] # mac
        self.splay_edge = [179, 182, 183]  # linux

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc."""
        materials = (
            MaterialGroup(tag=1, entities=self.slab_volume),
            MaterialGroup(tag=2, entities=self.wedge_volume),
            MaterialGroup(tag=3, entities=self.mantle_volume),
            MaterialGroup(tag=4, entities=self.crust_volume),
        )

        for material in materials:
            material.create_physical_group()

        face_groups = (
            BoundaryGroup(
                name="boundary_south", tag=10, dim=2, entities=self.surface_south
            ),
            BoundaryGroup(
                name="boundary_east", tag=11, dim=2, entities=self.surface_east
            ),
            BoundaryGroup(
                name="boundary_north", tag=12, dim=2, entities=self.surface_north
            ),
            BoundaryGroup(
                name="boundary_west", tag=13, dim=2, entities=self.surface_west
            ),
            BoundaryGroup(
                name="boundary_bottom", tag=14, dim=2, entities=self.surface_bottom
            ),
            BoundaryGroup(
                name="boundary_top", tag=15, dim=2, entities=self.surface_top
            ),
            BoundaryGroup(
                name="boundary_west_noslab",
                tag=16,
                dim=2,
                entities=self.surface_west_no_slab,
            ),
            BoundaryGroup(
                name="boundary_south_noslab",
                tag=17,
                dim=2,
                entities=self.surface_south_no_slab,
            ),
            BoundaryGroup(
                name="boundary_north_noslab",
                tag=18,
                dim=2,
                entities=self.surface_north_no_slab,
            ),
            BoundaryGroup(
                name="boundary_south_noslab_nosplay",
                tag=19,
                dim=2,
                entities=self.surface_south_no_slab_splay,
            ),
            BoundaryGroup(
                name="boundary_north_noslab_nosplay",
                tag=20,
                dim=2,
                entities=self.surface_north_no_slab_splay,
            ),
            BoundaryGroup(name="fault_slabtop", tag=30, dim=2, entities=self.slab_top),
            BoundaryGroup(
                name="fault_slabtop_edge", tag=130, dim=1, entities=self.slab_top_edge
            ),
            BoundaryGroup(
                name="fault_slabbot", tag=31, dim=2, entities=self.slab_bottom
            ),
            BoundaryGroup(
                name="fault_slabbot_edge",
                tag=131,
                dim=1,
                entities=self.slab_bottom_edge,
            ),
            BoundaryGroup(
                name="fault_slabtop_patch", tag=32, dim=2, entities=self.slab_top_patch
            ),
            BoundaryGroup(
                name="fault_slabtop_patch_edge",
                tag=132,
                dim=1,
                entities=self.slab_top_patch_edge,
            ),
            BoundaryGroup(
                name="fault_splay", tag=33, dim=2, entities=self.surface_splay
            ),
            BoundaryGroup(
                name="fault_splay_edge", tag=133, dim=1, entities=self.splay_edge
            ),
        )
        for group in face_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh."""
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # We setup a field `field_distance` with the distance from the top of the slab.
        field_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_distance, "SurfacesList", self.slab_top)

        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(
            field_distance, min_dx=self.DX_MIN, bias=self.DX_BIAS
        )
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        gmsh.model.mesh.generate(3)
        self.improve_quality()

    def _create_topography_surface(self):
        """Create topography/bathymetry surface."""
        dem = netCDF4.Dataset(self.FILENAME_LOCALDEM)
        latitude_1d = dem.variables["lat"][:]
        longitude_1d = dem.variables["lon"][:]
        elevation_2d = dem.variables["data"][:].astype(float)

        wires = []
        n_longitude = len(longitude_1d)
        n_latitude = len(latitude_1d)
        for i in range(n_latitude):
            points = numpy.zeros((n_longitude, 3), dtype=numpy.float64)
            points[:, 0] = latitude_1d[i]
            points[:, 1] = longitude_1d
            points[:, 2] = elevation_2d[i, :]

            coordsys.geoToMesh(points)
            wire = self._create_wire(points)
            wires.append(wire)
        topography_surface = gmsh.model.occ.add_thru_sections(
            wires, makeSolid=False, makeRuled=False, maxDegree=2
        )
        gmsh.model.occ.remove([(1, wire) for wire in wires], recursive=True)
        return topography_surface

    def _create_slab_top_surface(self):
        """Create surface for top of slab."""
        # Read slab contours
        contours = self._create_slab_splines()
        for contour in contours.values():
            coordsys.geoToMesh(contour)

        # Extend the slab to the topography
        contours_up_dip = self._create_extended_contours(contours)
        contours_all = dict(
            sorted((contours_up_dip | contours).items())
        )  # merge dicts and sort

        # Generate the top surface of the slab
        all_points_on_grid = []
        south_side = []
        north_side = []
        for contour in contours_all.values():
            gmsh_points = []
            for point in contour:
                gmsh_point = gmsh.model.occ.add_point(point[0], point[1], point[2])
                gmsh_points.append(gmsh_point)
            south_side.append(gmsh_points[0])
            north_side.append(gmsh_points[-1])
            all_points_on_grid.append(gmsh_points)

        spline_south = gmsh.model.occ.add_spline(south_side)
        spline_north = gmsh.model.occ.add_spline(north_side)
        spline_east = gmsh.model.occ.add_spline(all_points_on_grid[0])
        spline_west = gmsh.model.occ.add_spline(all_points_on_grid[-1])

        wire = gmsh.model.occ.add_wire(
            [
                spline_south,
                spline_east,
                spline_north,
                spline_west,
            ]
        )

        point_tags = numpy.concatenate(all_points_on_grid).tolist()
        surface = gmsh.model.occ.addSurfaceFilling(wire, pointTags=point_tags)
        gmsh.model.occ.remove([(0, point) for point in point_tags])
        return surface, contours

    def _create_slab_splines(self):
        """Read geographic coordinates from "cas_contours_dep.in.txt.gz" and return a dictionary where the keys are
        the depth of the contour and the values are the coordinates (latitude, longitude, elevation) as a numpy array.
        """
        with gzip.open(self.SLAB_CONTOURS_FILENAME, "rb") as file:
            lines = file.readlines()
        contours = {}
        points = []
        key = None
        for line in lines:
            if line.decode().strip() == "END":
                contours[key] = numpy.array(points, dtype=numpy.float64)
                points = []
            elif len(line.split()) == 1:
                key = int(line)
            else:
                pt = list(map(float, line.strip().split()))  # lon/lat/elev(km)
                this_point = (pt[1], pt[0], pt[2] * self.KM)  # lat/lon/elev(m)
                points.append(this_point)

        # The orientation of contour with key `5`` is reversed
        contours[5] = numpy.ascontiguousarray(numpy.flip(contours[5], axis=0))
        return contours

    def _calculate_local_strike(self, contour):
        """Calculate local strike along contour.
        """
        def _smooth_values(values):
            window_size = values.shape[0] // 3
            window = numpy.ones(window_size) / window_size

            pad_size = window_size // 2
            values_padded = numpy.pad(values, (pad_size, pad_size), mode="edge")
            values_smoothed = numpy.convolve(values_padded, window, mode="valid")
            return values_smoothed[0:-1]

        contour_smooth = numpy.stack(
            (_smooth_values(contour[:, 0]), _smooth_values(contour[:, 1]))
        ).T

        sx = contour_smooth[:-1, 0] - contour_smooth[1:, 0]
        sy = contour_smooth[:-1, 1] - contour_smooth[1:, 1]
        strike = numpy.atan2(sx, sy)
        strike = numpy.concatenate((numpy.array([strike[0]]), strike))
        return strike

    def _calculate_surface_avg_normal(self, surface):
        # Find the normals of the slab surface
        space = numpy.linspace(0, 1, 10)
        xv, yv = numpy.meshgrid(space, space)
        parametricCoords = numpy.column_stack((xv.flatten(), yv.flatten())).ravel()
        normals = gmsh.model.getNormal(surface, parametricCoords)
        normal = numpy.average(normals.reshape(-1, 3), axis=0)
        return normal

    def _create_extended_contours(self, contours):
        """Add contours up-dip from original contours.

        We increase the horizontal distance between the contours at a
        geometric rate. The first contour is at a distance of
        dist_horiz, followed by 2*dist_horiz, 4*dist_horiz, etc.

        The horizontal distance of contour n from the original one is
        (2**(n+1)-1)) * dist_horiz, n=0,1,2,...
        """
        n_contours = 5

        key = min(contours.keys())
        contour_top = contours[key]
        z_top = contour_top[0][2]
        dist_horiz = (self.UP_DIP_ELEV - z_top) / math.tan(self.UP_DIP_ANGLE)

        strike = self._calculate_local_strike(contour_top)
        dx = -dist_horiz * numpy.cos(strike)
        dy = dist_horiz * numpy.sin(strike)

        contours_up_dip = {}
        z_range = numpy.linspace(z_top, self.UP_DIP_ELEV, n_contours + 1)
        for i in range(n_contours):
            contour = numpy.array(contour_top)
            contour[:, 0] += (i + 1) * dx
            contour[:, 1] += (i + 1) * dy
            contour[:, 2] = z_range[i + 1]
            contours_up_dip[int(-z_range[i + 1] / self.KM)] = contour
        return contours_up_dip

    def _create_slab_volume(self, top_surface, surface_normal):
        """Create slab volume."""
        dx = self.SLAB_THICKNESS * surface_normal[0]
        dy = self.SLAB_THICKNESS * surface_normal[1]
        dz = self.SLAB_THICKNESS * surface_normal[2]
        tags = gmsh.model.occ.extrude(
            [(2, top_surface)], dx=dx, dy=dy, dz=dz, recombine=True
        )
        slab_volume = tags[1]
        return slab_volume

    def _create_splay_fault_surface(self, contours):
        """Create splay fault surface."""
        # Generate the splay fault
        splay_bottom = numpy.copy(contours[15])
        splay_top = numpy.copy(splay_bottom)

        splay_bottom[:, 2] -= 8.0e3
        splay_bottom_wire = self._create_wire(splay_bottom)

        splay_top[:, 0] -= 24.0e3
        splay_top[:, 2] = 3.0e3
        splay_top_wire = self._create_wire(splay_top)
        surface = gmsh.model.occ.add_thru_sections(
            [splay_bottom_wire, splay_top_wire],
            makeSolid=False,
            makeRuled=False,
            maxDegree=2,
        )
        return surface

    def _create_domain_box(self):
        """Create domain extending above topography/bathymetry."""
        # Generate the bounding box
        domain = gmsh.model.occ.add_box(
            -self.DOMAIN_X / 2,
            -self.DOMAIN_Y / 2,
            -self.DOMAIN_Z,
            self.DOMAIN_X,
            self.DOMAIN_Y,
            self.TOP_BOX + self.DOMAIN_Z,
        )
        return domain

    def _create_continental_moho_surface(self):
        """Create bottom of continental crust."""
        surface = self._add_xy_plane_surface_at_point(
            -self.DOMAIN_X / 2 - 50 * self.KM,
            -self.DOMAIN_Y / 2 - 50 * self.KM,
            -1 * self.KM
            - self.CRUST_DEPTH,  # Crust side of the bounding box about is approximately -1km below 0
            self.DOMAIN_X + 100 * self.KM,
            self.DOMAIN_Y + 100 * self.KM,
        )
        return surface

    @staticmethod
    def _create_wire(contour: numpy.ndarray):
        """Helper function for creating a spline and wire from a numpy array of shape (x,3)"""
        points = []
        for point in contour:
            gmsh_point = gmsh.model.occ.add_point(point[0], point[1], point[2])
            points.append(gmsh_point)
        spline = gmsh.model.occ.add_spline(points)
        wire = gmsh.model.occ.add_wire([spline])
        gmsh.model.occ.remove([(0, point) for point in points])
        # gmsh.model.occ.remove([(1, spline)])
        return wire

    @staticmethod
    def _add_xy_plane_surface_at_point(x, y, z, x_extent, y_extent):
        """Helper function for creating an xy-plane surface"""
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


# If script is called from the command line, run the application.
if __name__ == "__main__":
    App().main()
