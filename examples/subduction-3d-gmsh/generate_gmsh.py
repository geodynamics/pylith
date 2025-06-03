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

        positions = {}
        for depth, contour in contours.items():
            position = self._get_cartesian(contour[:,0],contour[:,1],contour[:,2])
            positions[depth] = position
        rotational_center = self._get_center(positions[5])
        rotation_matrix = self._get_rotation_matrix_from_direction(rotational_center)
        center = self._get_center(self._apply_rotation_points(positions[5], rotation_matrix))

        final_points = {}

        for depth,position in positions.items():
            position = self._apply_rotation_points(position, rotation_matrix)
            position = self._apply_centering_points(position, center)
            final_points[depth] = position

        splines = []
        wires = []

        sorted_final_points = [v for k, v in sorted(final_points.items())]
        for final_point in sorted_final_points:
            points = []
            for point in final_point:
                gmsh_point = gmsh.model.occ.add_point(point[0],point[1], point[2])
                points.append(gmsh_point)
            spline = gmsh.model.occ.add_spline(points)
            splines.append(spline)
            wires.append(gmsh.model.occ.add_wire([spline]))

        surface = gmsh.model.occ.add_thru_sections(wires,makeSolid=False, makeRuled=False)

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
    def _get_cartesian(lat_deg, lon_deg, alt):
        lat = np.radians(lat_deg)
        lon = np.radians(lon_deg)
        rad = np.float64(6378137.0)
        # Radius of the Earth (in meters)
        f = np.float64(1.0 / 298.257223563)  # Flattening factor WGS84 Model
        cosLat = np.cos(lat)
        sinLat = np.sin(lat)
        FF = (1.0 - f) ** 2
        C = 1 / np.sqrt(cosLat ** 2 + FF * sinLat ** 2)
        S = C * FF

        x = (rad * C + alt) * cosLat * np.cos(lon)
        y = (rad * C + alt) * cosLat * np.sin(lon)
        z = (rad * S + alt) * sinLat

        return np.vstack((x, y, z)).T

    @staticmethod
    def _get_center(vertices):
        return np.average(vertices, axis=0)

    @staticmethod
    def _get_rotation_matrix_from_direction(direction):
        normal = direction / np.linalg.norm(direction)  # Normalize the normal vector

        # Define the target normal vector (Z-axis)
        target = np.array([0, 0, 1])

        # Compute the rotation axis (cross product of normal and target)
        axis = np.cross(normal, target)
        axis = axis / np.linalg.norm(axis)  # Normalize the rotation axis

        # Compute the rotation angle (dot product of normal and target)
        angle = np.arccos(np.dot(normal, target))

        # Create the rotation matrix using Rodrigues' rotation formula
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
        R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)
        return R

    @staticmethod
    def _apply_rotation_points(points, R):
        return np.dot(points, R.T)

    @staticmethod
    def _apply_centering_points(points, center):
        return points - center

# If script is called from the command line, run the application.
if __name__ == "__main__":
    App().main()

# End of file

