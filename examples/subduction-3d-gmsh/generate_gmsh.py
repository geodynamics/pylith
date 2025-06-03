from operator import countOf

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

    UP_DIP_ELEV = 1.0*1000
    UP_DIP_DIST = 600.0*1000
    UP_DIP_ANGLE = 10.0*1000
    FAULT_STRIKE = 0.0
    CONTOURS_STRIDE = 4
    POINTS_STRIDE = 20

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
            points.append([pt[1], pt[0], pt[2]])  # lat/lon/elev
            all_points.append([pt[1], pt[0], pt[2]])
        return contours,np.array(all_points)

    def _generate_extended_contours(self,contours):
        key = min(contours.keys())
        contour_top = contours[key]
        z_top = contour_top[0][2] * 1000
        dist_horiz = (self.UP_DIP_ELEV - z_top) / math.tan(self.UP_DIP_ANGLE)
        dx = -dist_horiz * math.cos(self.FAULT_STRIKE)
        dy = dist_horiz * math.sin(self.FAULT_STRIKE)
        contours_up_dip = {}
        ncontours = int(math.ceil(math.log((self.UP_DIP_DIST / dist_horiz) + 1) / math.log(2.0)))
        for i in range(ncontours):
            contour = np.array(contour_top)
            contour[:, 0] += (2 ** i) * dx
            contour[:, 1] += (2 ** i) * dy
            contour[:, 2] = self.UP_DIP_ELEV
            contours_up_dip[-i] = contour
        return contours_up_dip

    def create_geometry(self):
        """Create geometry.
        """
        crs_wgs84_3d = CRS.from_epsg(4979)
        crs_projected = CRS.from_proj4(
            "+proj=tmerc +datum=WGS84 +lat_0=45.5231 +lon_0=-122.6765 +k=0.9996 +units=m +type=crs"
        )
        transformer = Transformer.from_crs(crs_wgs84_3d, crs_projected, always_xy=True)

        # gmsh.model.occ.add_box(-60*1000,-60*1000,-400*1000,800*1000,800*1000,400*1000)
        dem_nc = netCDF4.Dataset(self.FILENAME_LOCALDEM)
        latitude = dem_nc.variables["lat"][:]
        longitude = dem_nc.variables["lon"][:]
        elevation = dem_nc.variables["Band1"][:].astype(float)

        topography_points = np.full(elevation.shape,-1, dtype=np.int32)
        topo_mask = elevation.recordmask
        for i in range(elevation.shape[0]):
            for j in range(elevation.shape[1]):
                if topo_mask and topo_mask[i, j]:
                    continue  # Skip masked values
                point_latitude = latitude[i]
                point_longitude = longitude[j]
                point_elevation = elevation[i, j]
                position = transformer.transform(point_longitude, point_latitude,point_elevation)
                point = gmsh.model.occ.add_point(position[0], position[1], position[2])
                topography_points[i,j] = point


        contours,all_points = self._read_and_generate_splines()
        bounding_box = self._calculate_bounding_box(lats=all_points[:, 0], longs=all_points[:, 1])
        print("Bounding box:", bounding_box)
        contours[5] = np.flip(contours[5],axis=0)

        projected_contours = {}
        for depth, contour in contours.items():
            position = transformer.transform(contour[:,1],contour[:,0],contour[:,2]*1000)
            projected_contours[depth] = np.array(position).T

        # projected_contours_up_dip = self._generate_extended_contours(projected_contours)
        # projected_contours = projected_contours_up_dip | projected_contours

        splines = []
        wires = []

        sorted_project_contours = [v for k, v in sorted(projected_contours.items())]
        for projected_contour in sorted_project_contours:
            points = []
            for point in projected_contour:
                gmsh_point = gmsh.model.occ.add_point(point[0],point[1], point[2])
                points.append(gmsh_point)
            spline = gmsh.model.occ.add_spline(points)
            splines.append(spline)
            wires.append(gmsh.model.occ.add_wire([spline]))

        slab_top_surface = gmsh.model.occ.add_thru_sections(wires,makeSolid=False, makeRuled=False)
        gmsh.model.occ.synchronize()

        space = np.linspace(0, 1, 10)
        xv, yv = np.meshgrid(space, space)
        parametricCoords = np.column_stack((xv.flatten(), yv.flatten())).ravel()
        normals = gmsh.model.getNormal(slab_top_surface[0][1],parametricCoords)
        normal = np.average(normals.reshape(-1, 3),axis=0)

        dx = -self.SLAB_THICKNESS * normal[0]
        dy = -self.SLAB_THICKNESS * normal[1]
        dz = -self.SLAB_THICKNESS * normal[2]

        slab_volume = gmsh.model.occ.extrude(slab_top_surface,dx=dx,dy=dy,dz=dz,recombine=True)
        gmsh.model.occ.synchronize()
        gmsh.fltk.run()

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        pass

    def generate_mesh(self, cell):
        """Generate the mesh.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Laplace2D")

    @staticmethod
    def _calculate_bounding_box(lats, longs, radius=0):
        # Calculate the min and max for latitude and longitude
        min_lat = np.min(lats) - radius
        max_lat = np.max(lats) + radius
        min_lon = np.min(longs) - radius
        max_lon = np.max(longs) + radius

        return min_lat, max_lat, min_lon, max_lon

# If script is called from the command line, run the application.
if __name__ == "__main__":
    App().main()

# End of file

