from operator import countOf

import numpy as np
import gzip
# Import Gmsh Python interface
import gmsh

# Import the gmsh_utils Python module supplied with PyLith.
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)
from pyproj import Transformer, CRS

class App(GenerateMesh):
    """
    Application used to generate the mesh using Gmsh.

    App uses `GenerateMesh` from `gmsh_utils` for common functionality that we avoid
    duplicating in each of our examples.
    """
    CONTOURS_FILENAME = "cas_contours_dep.in.txt.gz"
    SLAB_THICKNESS = 50.0 * 1000
    SLAB_NORMAL_DIR = (+0.209, -0.016, +0.979)

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

    # def _create_points_from_file(self, filename):
    #     coordinates = numpy.loadtxt(filename)
    #     points = []
    #     for xy in coordinates:
    #         points.append(gmsh.model.occ.add_point(xy[0], xy[1], 0))
    #     return points

    def _read_and_generate_splines(self):
        with gzip.open(self.CONTOURS_FILENAME, "rb") as file:
            lines = file.readlines()
        contours = {}
        points = []
        key = None
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
        return contours

    def create_geometry(self):
        """Create geometry.
        """
        contours = self._read_and_generate_splines()
        contours[5] = np.flip(contours[5],axis=0)

        crs_wgs84_3d = CRS.from_epsg(4979)
        crs_projected = CRS.from_proj4(
            "+proj=tmerc +datum=WGS84 +lat_0=45.5231 +lon_0=-122.6765 +k=0.9996 +units=m +type=crs"
        )
        transformer = Transformer.from_crs(crs_wgs84_3d, crs_projected, always_xy=True)

        projected_contours = {}
        for depth, contour in contours.items():
            position = transformer.transform(contour[:,1],contour[:,0],contour[:,2]*1000)
            projected_contours[depth] = np.array(position).T

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

        dx = self.SLAB_THICKNESS * normal[0]
        dy = self.SLAB_THICKNESS * normal[1]
        dz = self.SLAB_THICKNESS * normal[2]

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

# If script is called from the command line, run the application.
if __name__ == "__main__":
    App().main()

# End of file

