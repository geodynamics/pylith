#!python3
"""WARNING: This script only works with the Python supplied with Cubit.

To run this script outside the Cubit GUI, you will need to run it like:

PATH_TO_CUBIT/Cubit.app/Contents/Frameworks/Python.framework/Versions/Current/bin/python3 generate_cubit.py --cubit-directory=PATH_TO_CUBIT/Cubit.app/Contents/MacOS

where you replace 'PATH_TO_CUBIT' with the absolute path.
"""


# -------------------------------------------------------------------------------------------------
# Utility functions. No edits should be needed.
# -------------------------------------------------------------------------------------------------
def setup_cubit():
    """Detect if we are running outside Cubit. If so, then we parse command line arguments and setup
    Cubit."""
    import sys

    if not "cubit" in sys.modules:
        import argparse
        import pathlib

        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--cubit-directory",
            action="store",
            dest="cubit_dir",
            required=True,
            help="Directory containing cubit executable.",
        )
        args = parser.parse_args()

        # Initialize cubit
        cubit_absdir = pathlib.Path(args.cubit_dir).expanduser().resolve()
        sys.path.append(str(cubit_absdir))
        import cubit

        cubit.init(["cubit", "-nojournal"])


def get_last_entity(entity_type: str):
    """Get last entity of type 'entity_type'.

    Example:
        domain = get_last_entity("volume")
    """
    get_type = {
        "vertex": cubit.vertex,
        "curve": cubit.curve,
        "surface": cubit.surface,
        "volume": cubit.volume,
        "body": cubit.body,
    }
    return get_type[entity_type](cubit.get_last_id(entity_type))


def get_entity(entity_type: str, name: str):
    """Get entity from name.

    Example:
        domain = get_entity("volume", name="domain")
    """
    get_type = {
        "vertex": cubit.vertex,
        "curve": cubit.curve,
        "surface": cubit.surface,
        "volume": cubit.volume,
        "body": cubit.body,
    }
    return get_type[entity_type](cubit.get_id_from_name(name))


def ids_from_entities(entities: list) -> str:
    """Create string with list of ids for 'entities'."""
    return cubit.get_id_string([entity.id() for entity in entities])


import gzip
import math
import sys
import os

import numpy

PYLITH_DIR = os.environ["PYLITH_DIR"]
sys.path.append(f"{PYLITH_DIR}/lib/python3.10/site-packages")

import coordsys

KM = 1000.0  # scale to convert km to m

SLAB_CONTOURS_FILENAME = "cas_contours_dep.in.txt.gz"
SURFACE_SLABTOP_FILENAME = "scratch/cubit_surf_splay.sat"
SLAB_THICKNESS = 50.0 * KM
DOMAIN_X = DOMAIN_Y = 800 * KM
DOMAIN_Z = 400 * KM

UP_DIP_ELEV = 5.0 * KM
UP_DIP_ANGLE = math.radians(30.0)

# Discretization size
DX_MIN = 20.0 * KM
DX_BIAS = 1.08
CELL = "tet"

# Depth of bottom of continental crust
CRUST_DEPTH = 30.0 * KM

# Size of patch on slab interface where we prescribe coseismic slip
PATCH_LENGTH = 200.0 * KM
PATCH_WIDTH = 275.0 * KM


# -------------------------------------------------------------------------------------------------
# Start cubit
# -------------------------------------------------------------------------------------------------
setup_cubit()
cubit.reset()


# -------------------------------------------------------------------------------------------------
# Geometry
# -------------------------------------------------------------------------------------------------
def _create_slab_points():
    """Read geographic coordinates from "cas_contours_dep.in.txt.gz" and return a dictionary where the keys are
    the depth of the contour and the values are the coordinates (latitude, longitude, elevation) as a numpy array.
    """
    with gzip.open(SLAB_CONTOURS_FILENAME, "rb") as file:
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
            this_point = (pt[1], pt[0], pt[2] * KM)  # lat/lon/elev(m)
            points.append(this_point)

    # The orientation of contour with key `5`` is reversed
    contours[5] = numpy.ascontiguousarray(numpy.flip(contours[5], axis=0))
    return contours


def _calculate_local_strike(contour):
    """Calculate local strike along contour."""

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
    strike = numpy.arctan2(sx, sy)
    strike = numpy.concatenate((numpy.array([strike[0]]), strike))
    return strike


def _calculate_surface_avg_normal(surface):
    """Calculate average surface normal."""
    # Find the normals of the slab surface
    space = numpy.linspace(0, 1, 10)
    xv, yv = numpy.meshgrid(space, space)
    parametricCoords = numpy.column_stack((xv.flatten(), yv.flatten()))
    normals = []
    for pointUV in parametricCoords:
        pointXYZ = surface.position_from_u_v(pointUV[0], pointUV[1])
        if pointXYZ[2] > 0.0:
            # Skip up-dip points
            continue
        normals.append(surface.normal_at(pointXYZ))
    normal = numpy.average(numpy.array(normals).reshape(-1, 3), axis=0)
    return normal


def _create_extended_contours(contours):
    """Add contours up-dip from original contours plus a few more westward."""
    n_contours_updip = 5
    n_contours_horiz = 2

    key = min(contours.keys())
    contour_top = contours[key]
    z_top = contour_top[0][2]
    dist_horiz = (UP_DIP_ELEV - z_top) / math.tan(UP_DIP_ANGLE)

    strike = _calculate_local_strike(contour_top)
    dx = -dist_horiz * numpy.cos(strike)
    dy = dist_horiz * numpy.sin(strike)

    contours_up_dip = {}
    z_range = numpy.linspace(z_top, UP_DIP_ELEV, n_contours_updip + 1)
    for i in range(n_contours_updip):
        contour = numpy.array(contour_top)
        contour[:, 0] += (i + 1) * dx
        contour[:, 1] += (i + 1) * dy
        contour[:, 2] = z_range[i + 1]
        contours_up_dip[int(-z_range[i + 1] / KM)] = contour
    for i in range(n_contours_horiz):
        contour = numpy.array(contours_up_dip[int(-UP_DIP_ELEV / KM)])
        contour[:, 0] += -80.0 * (i + 1) * KM
        contours_up_dip[int(-(i + 2) * UP_DIP_ELEV / KM)] = contour
    return contours_up_dip


def _create_slab_top_surface():
    """Create top surface of slab."""
    # Read slab contours
    contours = _create_slab_points()
    for contour in contours.values():
        coordsys.geoToMesh(contour)

    # Extend the slab to the topography
    contours_up_dip = _create_extended_contours(contours)
    contours_all = dict(
        sorted(list(contours_up_dip.items()) + list(contours.items()))
    )  # merge dicts and sort

    # Generate the top surface of the slab
    splines = []
    for contour in contours_all.values():
        cubit_points = []
        for point in contour:
            cubit_point = cubit.create_vertex(point[0], point[1], point[2])
            cubit_points.append(cubit_point)
        cubit.cmd(f"create curve spline vertex {ids_from_entities(cubit_points)}")
        spline = get_last_entity("curve")
        splines.append(spline)

    cubit.cmd(f"create surface skin curve {ids_from_entities(splines)}")
    cubit.cmd("delete curve all")
    cubit.cmd("delete vertex all")
    cubit.cmd(f"export acis '{SURFACE_SLABTOP_FILENAME}' overwrite")
    surface = get_last_entity("surface")
    return surface, contours


def _create_slab_volume(top_surface, surface_normal):
    """Create slab volume."""
    dx = SLAB_THICKNESS * -surface_normal[0]
    dy = SLAB_THICKNESS * -surface_normal[1]
    dz = SLAB_THICKNESS * -surface_normal[2]
    cubit.cmd(f"surface {top_surface.id()} move direction {dx} {dy} {dz} copy")
    s_slabbot = get_last_entity("surface")
    cubit.cmd(f"create body loft surface {top_surface.id()} {s_slabbot.id()}")
    slab_volume = get_last_entity("volume")
    cubit.cmd(f"delete surface {s_slabbot.id()}")
    return slab_volume


def _create_splay_fault_surface(contours):
    """Create splay fault surface from slab contours.

    We use one slab contour for the intersection of the splay fault with the top of the slab.
    """
    # Generate the splay fault
    splay_bottom = numpy.copy(contours[15])
    splay_top = numpy.copy(splay_bottom)

    splay_bottom[:, 2] -= 8.0e3

    splay_top[:, 0] -= 24.0e3
    splay_top[:, 2] = 3.0e3

    splines = []
    for contour in [splay_bottom, splay_top]:
        cubit_points = []
        for point in contour:
            cubit_point = cubit.create_vertex(point[0], point[1], point[2])
            cubit_points.append(cubit_point)
        cubit.cmd(f"create curve spline vertex {ids_from_entities(cubit_points)}")
        spline = get_last_entity("curve")
        splines.append(spline)

    cubit.cmd(f"create surface skin curve {ids_from_entities(splines)}")
    cubit.cmd("delete curve all")
    cubit.cmd("delete vertex all")
    cubit.cmd(f"export acis '{SURFACE_SLABTOP_FILENAME}' overwrite")
    surface = get_last_entity("surface")
    return surface


# Create domain
v_domain = cubit.brick(DOMAIN_X, DOMAIN_Y, DOMAIN_Z)
v_domain.remove_entity_names()
v_domain.set_entity_name("domain")
cubit.move(v_domain, (-25 * KM, 0.0, -0.5 * DOMAIN_Z))

# Create top surface of slab.
s_slabtop, contours = _create_slab_top_surface()
slab_normal = _calculate_surface_avg_normal(s_slabtop)

# Create slab volume.
v_slab = _create_slab_volume(s_slabtop, slab_normal)
v_slab.remove_entity_names()
v_slab.set_entity_name("slab")

# Create splay fault.
s_splay_fault = _create_splay_fault_surface(contours)
s_splay_fault.remove_entity_names()
s_splay_fault.set_entity_name("splay_fault")

# Create patch for top of slab interface.
v_patch = cubit.brick(PATCH_WIDTH, PATCH_LENGTH, 200 * KM)
v_patch.remove_entity_names()
v_patch.set_entity_name("patch")
cubit.move(v_patch, (-DOMAIN_X / 2 + PATCH_WIDTH / 2, 0, 0))

# Create continental moho (bottom of continental crust)
cubit.cmd(f"create planar surface with plane zplane offset {-CRUST_DEPTH}")
s_moho = get_last_entity("surface")
s_moho.remove_entity_names()
s_moho.set_entity_name("moho")

# Chop domain with slab
cubit.cmd(f"chop volume {v_domain.id()} with volume {v_slab.id()}")
v_slab = get_entity("volume", "domain")
v_slab.remove_entity_names()
v_slab.set_entity_name("slab")
v_domain_other = get_entity("volume", "domain@A")

# Chop domain (continent) with moho
cubit.cmd(f"webcut volume {v_domain_other.id()} with plane surface {s_moho.id()}")
cubit.cmd(f"delete surface {s_moho.id()}")
v_crust = get_entity("volume", "domain@A")
v_crust.remove_entity_names()
v_crust.set_entity_name("crust")
v_mantle = get_entity("volume", "domain@B")
v_mantle.remove_entity_names()
v_mantle.set_entity_name("mantle")

# Chop crust with splay fault
cubit.cmd(f"webcut volume {v_crust.id()} with sheet surface {s_splay_fault.id()}")
cubit.cmd(f"delete body {ids_from_entities(s_splay_fault.bodies())}")
v_wedge = get_entity("volume", "crust")
v_wedge.remove_entity_names()
v_wedge.set_entity_name("wedge")
v_crust = get_entity("volume", "crust@A")
v_crust.remove_entity_names()
v_crust.set_entity_name("crust")

# Enscribe the rupture patch onto the geometry of the subducting slab
# surface by chopping the original slab volume and then imprinting the surface
# onto the wedge and crust.
v_slabtop = s_slabtop.volumes()[0]
v_slabtop.remove_entity_names()
v_slabtop.set_entity_name("slabtop")
cubit.cmd(f"chop volume {v_slabtop.id()} with volume {v_patch.id()}")
v_slab_patch = get_entity("volume", "slabtop")
v_slab_other = get_entity("volume", "slabtop@A")
cubit.cmd(f"imprint volume {v_wedge.id()} with volume {v_slab_patch.id()}")
cubit.cmd(f"imprint volume {v_crust.id()} with volume {v_slab_patch.id()}")
cubit.cmd(f"delete volume {v_slab_patch.id()}")
cubit.cmd(f"delete volume {v_slab_other.id()}")


cubit.cmd("imprint all")
cubit.cmd("merge all")

# Manually get surface ids from GUI
# idless = cubit.get_idless_signatures("surface", [29, 43, 49, 56])
s_north_str = "Surface ( at -171362 400000 -73929.7 ordinal 1 )  Surface ( at -25000 400000 -215000 ordinal 1 )  Surface ( at -347063 400000 -6763 ordinal 1 )  Surface ( at 32505.5 400000 -15000 ordinal 1 )"
s_north_noslab_str = "Surface ( at -25000 400000 -215000 ordinal 1 )  Surface ( at -347063 400000 -6763 ordinal 1 )  Surface ( at 32505.5 400000 -15000 ordinal 1 )"
s_north_noslab_nosplay_str = "Surface ( at -25000 400000 -215000 ordinal 1 )  Surface ( at 32505.5 400000 -15000 ordinal 1 )"
# idless = cubit.get_idless_signatures("surface", [28, 36])
s_west_str = "Surface ( at -425000 8.73115e-11 -23145.1 ordinal 1 )  Surface ( at -425000 8.73115e-11 -220874 ordinal 1 )"
s_west_noslab_str = "Surface ( at -425000 8.73115e-11 -220874 ordinal 1 )"
# idless = cubit.get_idless_signatures("surface", [26, 46, 51, 54])
s_south_str = "Surface ( at -140909 -400000 -73929.7 ordinal 1 )  Surface ( at -25000 -400000 -215000 ordinal 1 )  Surface ( at -186041 -400000 -6439.4 ordinal 1 )  Surface ( at 114076 -400000 -15000 ordinal 1 ) "
s_south_noslab_str = "Surface ( at -25000 -400000 -215000 ordinal 1 )  Surface ( at -186041 -400000 -6439.4 ordinal 1 )  Surface ( at 114076 -400000 -15000 ordinal 1 ) "
s_south_noslab_nosplay_str = "Surface ( at -25000 -400000 -215000 ordinal 1 )  Surface ( at 114076 -400000 -15000 ordinal 1 ) "
# idless = cubit.get_idless_signatures("surface", [40, 44])
s_east_str = "Surface ( at 375000 0 -15000 ordinal 1 )  Surface ( at 375000 0 -215000 ordinal 1 )"
# idless = cubit.get_idless_signatures("surface", [2])
s_bottom_str = "Surface ( at -25000 0 -400000 ordinal 1 )"
# idless = cubit.get_idless_signatures("surface", [27, 48, 53])
s_top_str = "Surface ( at -331061 -4.65661e-10 0 ordinal 1 )  Surface ( at -271736 -2.91038e-11 0 ordinal 1 )  Surface ( at 32505.5 -2.91038e-11 0 ordinal 1 )"

# idless = cubit.get_idless_signatures("surface", [45, 59, 60, 61, 65, 66])
s_fault_slabtop_str = "Surface ( at -19644.1 35477.7 -43492.6 ordinal 1 )  Surface ( at -262020 10520.6 -998.027 ordinal 1 )  Surface ( at -262020 10520.6 -998.027 ordinal 2 )  Surface ( at -262020 10520.6 -998.027 ordinal 3 )  Surface ( at -178716 -14836.1 -11135 ordinal 1 )  Surface ( at -178716 -14836.1 -11135 ordinal 2 )"
# idless = cubit.get_idless_signatures("surface", [23])
s_fault_slabbot_str = "Surface ( at -218964 43022.7 -55511.1 ordinal 1 )"
# idless = cubit.get_idless_signatures("surface", [47])
s_fault_splay_str = "Surface ( at -178906 -634.119 3000 ordinal 1 )"
# idless = cubit.get_idless_signatures("surface", [61, 65])
s_fault_slabtop_patch_str = "Surface ( at -262020 10520.6 -998.027 ordinal 3 )  Surface ( at -178716 -14836.1 -11135 ordinal 1 )"

# Create lists of curve ids for fault buried edges
# idless = cubit.get_idless_signatures("curve", [87])
c_fault_slabtop_str = "Curve ( at 141358 5406.16 -100000 ordinal 1 )  "
# idless = cubit.get_idless_signatures("curve", [85])
c_fault_slabbot_str = "Curve ( at 129041 5568.68 -147859 ordinal 1 )"
# idless = cubit.get_idless_signatures("curve", [153, 155, 157])
c_fault_splay_str = "Curve ( at -202255 262739 -13594.3 ordinal 1 )  Curve ( at -155213 -250596 -13452.1 ordinal 1 )  Curve ( at -163629 21.9472 -13539.2 ordinal 1 )"
# idless = cubit.get_idless_signatures("curve", [154, 169, 168, 167, 151])
c_fault_slabtop_patch_str = "Curve ( at -214293 -100000 -6071.25 ordinal 1 )  Curve ( at -143578 -100000 -16859.5 ordinal 1 )  Curve ( at -125000 -0.0682128 -20910 ordinal 1 )  Curve ( at -145033 100000 -16952.5 ordinal 1 )  Curve ( at -222375 100000 -6445.98 ordinal 1 )"


# -------------------------------------------------------------------------------------------------
# Generate the mesh
# -------------------------------------------------------------------------------------------------
cubit.cmd("volume all scheme tetmesh")

# Reset sizes
cubit.cmd("curve all scheme default")
cubit.cmd("surface all sizing function none")
cubit.cmd("volume all sizing function none")

# Set size on faults
cubit.cmd(f"surface {s_fault_slabtop_str} size {DX_MIN}")

cubit.cmd(
    f"volume all sizing function skeleton min_size {DX_MIN} max_gradient {DX_BIAS} max_size 50e+3"
)

# cubit.cmd("preview mesh volume all")
cubit.cmd("set tetmesher optimize sliver on")
cubit.cmd("mesh volume all")

# Smooth the mesh to improve mesh quality
# :NOTE: We have commented out smoothing the mesh, because there is a bug in Cubit v16.16 which
# causes Cubit to crash for some commands after smoothing.
# cubit.cmd("volume all smooth scheme condition number beta 1.7 cpu 10")
# cubit.cmd("smooth volume all")


# -------------------------------------------------------------------------------------------------
# Create blocks for materials and nodesets for boundary conditions.
# -------------------------------------------------------------------------------------------------
# We follow the general approach of creating groups and then creating the nodesets
# from the groups, so that we can apply boolean operations (e.g., union, intersection)
# on the groups to create the desired nodesets.

# Blocks
cubit.cmd(f"block 1 volume {v_slab.id()}")
cubit.cmd("block 1 name 'slab'")

cubit.cmd(f"block 2 volume {v_wedge.id()}")
cubit.cmd("block 2 name 'wedge'")

cubit.cmd(f"block 3 volume {v_mantle.id()}")
cubit.cmd("block 3 name 'mantle'")

cubit.cmd(f"block 4 volume {v_crust.id()}")
cubit.cmd("block 4 name 'crust'")


# Starting in PyLith v5, we can use sidesets instead of nodesets for BCs.

# Create sideset for south boundary
cubit.cmd(f"group 'boundary_south' add {s_south_str}")
cubit.cmd("sideset 10 group boundary_south")
cubit.cmd("sideset 10 name 'boundary_south'")

# Create sideset for east boundary
cubit.cmd(f"group 'boundary_east' add {s_east_str}")
cubit.cmd("sideset 11 group boundary_east")
cubit.cmd("sideset 11 name 'boundary_east'")

# Create sideset for north boundary
cubit.cmd(f"group 'boundary_north' add {s_north_str}")
cubit.cmd("sideset 12 group boundary_north")
cubit.cmd("sideset 12 name 'boundary_north'")

# Create sideset for west boundary
cubit.cmd(f"group 'boundary_west' add {s_west_str}")
cubit.cmd("sideset 13 group boundary_west")
cubit.cmd("sideset 13 name 'boundary_west'")

# Create sideset for bottom boundary
cubit.cmd(f"group 'boundary_bottom' add {s_bottom_str}")
cubit.cmd("sideset 14 group boundary_bottom")
cubit.cmd("sideset 14 name 'boundary_bottom'")

# Create sideset for top boundary
cubit.cmd(f"group 'boundary_top' add {s_top_str}")
cubit.cmd("sideset 15 group boundary_top")
cubit.cmd("sideset 15 name 'boundary_top'")


# Boundaries without slab
cubit.cmd(f"group 'boundary_west_noslab' add {s_west_noslab_str}")
cubit.cmd("sideset 16 group boundary_west_noslab")
cubit.cmd("sideset 16 name 'boundary_west_noslab'")

cubit.cmd(f"group 'boundary_south_noslab' add {s_south_noslab_str}")
cubit.cmd("sideset 17 group boundary_south_noslab")
cubit.cmd("sideset 17 name 'boundary_south_noslab'")

cubit.cmd(f"group 'boundary_north_noslab' add {s_north_noslab_str}")
cubit.cmd("sideset 18 group boundary_north_noslab")
cubit.cmd("sideset 18 name 'boundary_north_noslab'")

# Boundaries without slab and without splay
cubit.cmd(f"group 'boundary_south_noslab_nosplay' add {s_south_noslab_nosplay_str}")
cubit.cmd("sideset 19 group boundary_south_noslab_nosplay")
cubit.cmd("sideset 19 name 'boundary_south_noslab_nosplay'")

cubit.cmd(f"group 'boundary_north_noslab_nosplay' add {s_north_noslab_nosplay_str}")
cubit.cmd("sideset 20 group boundary_north_noslab_nosplay")
cubit.cmd("sideset 20 name 'boundary_north_noslab_nosplay'")


# Create sideset for fault 'slabtop'
cubit.cmd(f"group 'fault_slabtop' add {s_fault_slabtop_str}")
cubit.cmd("sideset 30 group fault_slabtop")
cubit.cmd("sideset 30 name 'fault_slabtop'")

# Create sideset for fault 'slabbot'
cubit.cmd(f"group 'fault_slabbot' add {s_fault_slabbot_str}")
cubit.cmd("sideset 31 group fault_slabbot")
cubit.cmd("sideset 31 name 'fault_slabbot'")

# Create sideset for fault 'slabtop_patch'
cubit.cmd(f"group 'fault_slabtop_patch' add {s_fault_slabtop_patch_str}")
cubit.cmd("sideset 32 group fault_slabtop_patch")
cubit.cmd("sideset 32 name 'fault_slabtop_patch'")

# Create sideset for fault 'splay'
cubit.cmd(f"group 'fault_splay' add {s_fault_splay_str}")
cubit.cmd("sideset 33 group fault_splay")
cubit.cmd("sideset 33 name 'fault_splay'")

# Buried edges (nodesets)

# Create nodeset for fault 'slabtop' edges
cubit.cmd(f"group 'fault_slabtop_edge' add node in {c_fault_slabtop_str}")
cubit.cmd("nodeset 130 group fault_slabtop_edge")
cubit.cmd("nodeset 130 name 'fault_slabtop_edge'")

# Create nodeset for fault 'slabbot' edges
cubit.cmd(f"group 'fault_slabbot_edge' add node in {c_fault_slabbot_str}")
cubit.cmd("nodeset 131 group fault_slabbot_edge")
cubit.cmd("nodeset 131 name 'fault_slabbot_edge'")

# Create nodeset for fault 'slabtop_patch' edges
cubit.cmd(f"group 'fault_slabtop_patch_edge' add node in {c_fault_slabtop_patch_str}")
cubit.cmd("nodeset 132 group fault_slabtop_patch_edge")
cubit.cmd("nodeset 132 name 'fault_slabtop_patch_edge'")

# Create nodeset for fault 'splay' edges
cubit.cmd(f"group 'fault_splay_edge' add node in {c_fault_splay_str}")
cubit.cmd("nodeset 133 group fault_splay_edge")
cubit.cmd("nodeset 133 name 'fault_splay_edge'")


# Write mesh as ExodusII file
cubit.cmd(f"export mesh 'input/mesh_{CELL}.exo' dimension 3 overwrite")
