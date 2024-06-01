#!python3
"""WARNING: This script only works with Python supplied with Cubit.

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
        parser.add_argument("--cubit-directory", action="store", dest="cubit_dir", required=True, help="Directory containing cubit executable.")
        args = parser.parse_args()

        # Initialize cubit
        cubit_absdir = pathlib.Path(args.cubit_dir).expanduser().resolve()
        sys.path.append(str(cubit_absdir))
        import cubit
        cubit.init(['cubit','-nojournal'])


def get_last_entity(entity_type: str):
    """Get last entry of type 'entity_type'.
    
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


def ids_from_entities(entities: list) -> str:
    """Get string with list of ids for 'entities'."""
    return cubit.get_id_string([entity.id() for entity in entities])


def create_points_from_file(filename: str) -> list:
    """Create points (vertices) from file with coordinates of points in mesh coordinate system."""
    coordinates = numpy.loadtxt(filename)
    points = []
    for xy in coordinates:
        points.append(cubit.create_vertex(xy[0], xy[1], 0))
    return points


# -------------------------------------------------------------------------------------------------
# Mesh definition
# -------------------------------------------------------------------------------------------------
#
#    pNW---------------------------------pNE
#    |                                   |
#    |                                   |
#    |         main                      |
#    |          \          east          |
#    |           \        /              |
#    |            \      /               |
#    |             \    /                |
#    |              pI                   |
#    |             /   \                 |
#    |            /     \                |
#    |          /        \               |
#    |        /                          |
#    |      /                            |
#    |     west                          |
#    |                                   |
#    |                                   |
#    pSW---------------------------------pSE
#
# The fault traces intersect at point pI.
#
# We create a rectangular prism for the domain and create fault surfaces from the surface traces.
# We then subdivide the domain using the fault surfaces. With this workflow we keep track of
# the ids of objects rather than naming them.

import numpy

# Locations defining domain
km = 1000.0
DOMAIN_CENTER = (453700.0, 3947000.0)
DOMAIN_X = 80*km
DOMAIN_Y = 60*km
DOMAIN_Z = 40*km
FAULT_DEPTH = 15.0*km

# Files with coordinates of fault traces in UTM zone 11
FILENAME_MAINTRACE = "faulttrace_main_utm.txt"
FILENAME_WESTTRACE = "faulttrace_west_utm.txt"
FILENAME_STRANDETRACE = "faulttrace_east_utm.txt"

# Discretization size on faults and bias (rate of cell size increase away from fault)
DX_FAULT = 2.5*km
DX_BIAS = 1.07

CELL = "tet"

# -------------------------------------------------------------------------------------------------
# Start cubit
# -------------------------------------------------------------------------------------------------
setup_cubit()
cubit.reset()

# -------------------------------------------------------------------------------------------------
# Geometry
# -------------------------------------------------------------------------------------------------
# Create domain
domain = cubit.brick(DOMAIN_X, DOMAIN_Y, DOMAIN_Z)
cubit.move(domain, (DOMAIN_CENTER[0], DOMAIN_CENTER[1], -0.5*DOMAIN_Z))


# Create main fault
points_fault_main = create_points_from_file(FILENAME_MAINTRACE)
cubit.cmd(f"create curve spline vertex {points_fault_main[0].id()} to {points_fault_main[-1].id()}")
faulttrace = get_last_entity("curve")
cubit.cmd(f"sweep curve {faulttrace.id()} vector 0 0 -1 distance {FAULT_DEPTH}")
surface_fault_main = get_last_entity("surface")

# Create points for west fault strand
points_fault_west = create_points_from_file(FILENAME_WESTTRACE)
cubit.cmd(f"create curve spline vertex {points_fault_west[0].id()} to {points_fault_west[-1].id()}")
faulttrace = get_last_entity("curve")
cubit.cmd(f"sweep curve {faulttrace.id()} vector 0 0 -1 distance {FAULT_DEPTH}")
surface_fault_west = get_last_entity("surface")

# Create points for east fault strand
points_fault_east = create_points_from_file(FILENAME_STRANDETRACE)
cubit.cmd(f"create curve spline vertex {points_fault_east[0].id()} to {points_fault_east[-1].id()}")
faulttrace = get_last_entity("curve")
cubit.cmd(f"sweep curve {faulttrace.id()} vector 0 0 -1 distance {FAULT_DEPTH}")
surface_fault_east = get_last_entity("surface")

# Remove fault trace vertices now that we have fault surfaces
cubit.cmd("delete vertex all")


# Create extended surface for fault 'main' and update ids for domain (split by fault)
cubit.cmd(f"webcut volume {domain.id()} with sheet extended from surface {surface_fault_main.id()}")
domain_west = get_last_entity("volume")
domain_east = domain

# Create extended surface for fault 'west' and update ids for domain (split by fault)
cubit.cmd(f"webcut volume {domain_west.id()} with sheet extended from surface {surface_fault_west.id()}")
domain_nw = domain_west
domain_sw = get_last_entity("volume")
del domain_west

# Create extended surface for fault 'east' and update ids for domain (split by fault)
cubit.cmd(f"webcut volume {domain_east.id()} with sheet extended from surface {surface_fault_east.id()}")
domain_ne = domain_east
domain_se = get_last_entity("volume")
del domain_east
domain_ids = ids_from_entities([domain_nw, domain_sw, domain_se, domain_ne])

# Imprint/merge
fault_curve_ids = ids_from_entities(surface_fault_main.curves() + surface_fault_west.curves() + surface_fault_east.curves())
cubit.cmd(f"imprint volume {domain_ids} with curve {fault_curve_ids}")
cubit.cmd(f"delete surface {surface_fault_main.id()} {surface_fault_west.id()} {surface_fault_east.id()}")
cubit.cmd("imprint all")
cubit.cmd("merge all")

# Manually get surface ids from GUI
#idless = cubit.get_idless_signatures("surface", [17, 31, 36])
surfaces_north_str = "Surface ( at 427180 3.977e+06 -20000 ordinal 1 )  Surface ( at 455807 3.977e+06 -20000 ordinal 1 )  Surface ( at 482327 3.977e+06 -20000 ordinal 1 )"
#idless = cubit.get_idless_signatures("surface", [4])
surfaces_west_str = "Surface ( at 413700 3.947e+06 -20000 ordinal 1 )"
#idless = cubit.get_idless_signatures("surface", [12, 21, 26])
surfaces_south_str = "Surface ( at 487831 3.917e+06 -20000 ordinal 1 )  Surface ( at 417569 3.917e+06 -20000 ordinal 1 )  Surface ( at 451700 3.917e+06 -20000 ordinal 1 )"
#idless = cubit.get_idless_signatures("surface", [6])
surfaces_east_str = "Surface ( at 493700 3.947e+06 -20000 ordinal 1 )"
#idless = cubit.get_idless_signatures("surface", [24, 27, 32, 39])
surfaces_bottom_str = "Surface ( at 433718 3.947e+06 -40000 ordinal 1 )  Surface ( at 451700 3.93223e+06 -40000 ordinal 1 )  Surface ( at 455807 3.96223e+06 -40000 ordinal 1 )  Surface ( at 473718 3.947e+06 -40000 ordinal 1 )"
#idless = cubit.get_idless_signatures("surface", [22, 29, 34, 37])
surfaces_top_str = "Surface ( at 433718 3.947e+06 0 ordinal 1 )  Surface ( at 451700 3.93223e+06 0 ordinal 1 )  Surface ( at 455807 3.96223e+06 0 ordinal 1 )  Surface ( at 473718 3.947e+06 0 ordinal 1 )"

#idless = cubit.get_idless_signatures("surface", [43, 51])
surfaces_fault_main_str = "Surface ( at 444628 3.96154e+06 -20000 ordinal 2 )  Surface ( at 468984 3.93314e+06 -20000 ordinal 2 )"
#idless = cubit.get_idless_signatures("surface", [46])
surfaces_fault_west_str = "Surface ( at 437634 3.93218e+06 -20000 ordinal 2 )"
#idless = cubit.get_idless_signatures("surface", [40])
surfaces_fault_east_str = "Surface ( at 462524 3.96213e+06 -20000 ordinal 1 )"

# Create lists of curve ids for fault buried edges
#idless = cubit.get_idless_signatures("curve", [77, 80, 97, 100])
fault_main_edges_str = "Curve ( at 444204 3.96317e+06 -7500 ordinal 1 )  Curve ( at 447417 3.95436e+06 -15000 ordinal 1 )  Curve ( at 466418 3.93633e+06 -7500 ordinal 1 )  Curve ( at 460642 3.94255e+06 -15000 ordinal 1 )"
#idless = cubit.get_idless_signatures("curve", [83, 86, 75])
fault_west_edges_str = "Curve ( at 445031 3.93911e+06 -7500 ordinal 1 )  Curve ( at 449604 3.94307e+06 -15000 ordinal 1 )  Curve ( at 453735 3.94746e+06 -7500 ordinal 1 )"
#idless = cubit.get_idless_signatures("curve", [73, 75, 76])
fault_east_edges_str = "Curve ( at 457036 3.95244e+06 -7500 ordinal 1 )  Curve ( at 453735 3.94746e+06 -7500 ordinal 1 )  Curve ( at 455387 3.94994e+06 -15000 ordinal 1 )"


# -------------------------------------------------------------------------------------------------
# Generate the mesh
# -------------------------------------------------------------------------------------------------
cubit.cmd("volume all scheme tetmesh")

# Reset sizes
cubit.cmd("curve all scheme default")
cubit.cmd("surface all sizing function none")
cubit.cmd("volume all sizing function none")

# Set size on faults
faults_str = f"{surfaces_fault_main_str} {surfaces_fault_west_str} {surfaces_fault_east_str}"
cubit.cmd(f"surface {faults_str} size {DX_FAULT}")

cubit.cmd(f"volume all sizing function skeleton min_size {DX_FAULT} max_gradient {DX_BIAS}")

# cubit.cmd("preview mesh volume all")
cubit.cmd("set tetmesher optimize sliver on")
cubit.cmd("mesh volume all")

# Smooth the mesh to improve mesh quality
# :NOTE: We have commented out smoothing the mesh, because there is a bug in Cubit v16.16 which
# causes Cubit to crash for some commands after smoothing.
#cubit.cmd("volume all smooth scheme condition number beta 1.7 cpu 10")
#cubit.cmd("smooth volume all")


# -------------------------------------------------------------------------------------------------
# Create blocks for materials and nodesets for boundary conditions.
# -------------------------------------------------------------------------------------------------
# We follow the general approach of creating groups and then creating the nodesets
# from the groups, so that we can apply boolean operations (e.g., union, intersection)
# on the groups to create the desired nodesets.

# Blocks
cubit.cmd(f"block 1 volume {domain_ids}")
cubit.cmd("block 1 name 'domain'")


# Nodesets

if True:
    # Create nodeset for south boundary
    cubit.cmd(f"group 'boundary_south' add node in {surfaces_south_str}")
    cubit.cmd("nodeset 10 group boundary_south")
    cubit.cmd("nodeset 10 name 'boundary_south'")

    # Create nodeset for east boundary
    cubit.cmd(f"group 'boundary_east' add node in {surfaces_east_str}")
    cubit.cmd("nodeset 11 group boundary_east")
    cubit.cmd("nodeset 11 name 'boundary_east'")

    # Create nodeset for north boundary
    cubit.cmd(f"group 'boundary_north' add node in {surfaces_north_str}")
    cubit.cmd("nodeset 12 group boundary_north")
    cubit.cmd("nodeset 12 name 'boundary_north'")

    # Create nodeset for west boundary
    cubit.cmd(f"group 'boundary_west' add node in {surfaces_west_str}")
    cubit.cmd("nodeset 13 group boundary_west")
    cubit.cmd("nodeset 13 name 'boundary_west'")

    # Create nodeset for bottom boundary
    cubit.cmd(f"group 'boundary_bottom' add node in {surfaces_bottom_str}")
    cubit.cmd("nodeset 14 group boundary_bottom")
    cubit.cmd("nodeset 14 name 'boundary_bottom'")

    # Create nodeset for top boundary
    cubit.cmd(f"group 'boundary_top' add node in {surfaces_top_str}")
    cubit.cmd("nodeset 15 group boundary_top")
    cubit.cmd("nodeset 15 name 'boundary_top'")


    # Create nodeset for fault 'main'
    cubit.cmd(f"group 'fault_main' add node in {surfaces_fault_main_str}")
    cubit.cmd("nodeset 20 group fault_main")
    cubit.cmd("nodeset 20 name 'fault_main'")

    # Create nodeset for fault 'main' edges
    cubit.cmd(f"group 'fault_main_edges' add node in {fault_main_edges_str}")
    cubit.cmd("nodeset 30 group fault_main_edges")
    cubit.cmd("nodeset 30 name 'fault_main_edges'")

    # Create nodeset for fault 'west'
    cubit.cmd(f"group 'fault_west' add node in {surfaces_fault_west_str}")
    cubit.cmd("nodeset 21 group fault_west")
    cubit.cmd("nodeset 21 name 'fault_west'")

    # Create nodeset for fault 'west' edges
    cubit.cmd(f"group 'fault_west_edges' add node in {fault_west_edges_str}")
    cubit.cmd("nodeset 31 group fault_west_edges")
    cubit.cmd("nodeset 31 name 'fault_west_edges'")

    # Create nodeset for fault 'east'
    cubit.cmd(f"group 'fault_east' add node in {surfaces_fault_east_str}")
    cubit.cmd("nodeset 22 group fault_east")
    cubit.cmd("nodeset 22 name 'fault_east'")

    # Create nodeset for fault 'east' edges
    cubit.cmd(f"group 'fault_east_edges' add node in {fault_east_edges_str}")
    cubit.cmd("nodeset 32 group fault_east_edges")
    cubit.cmd("nodeset 32 name 'fault_east_edges'")


# Starting in PyLith v5, we will use sidesets instead of nodesets for BCs.
# Sidesets
#
if False:
    # Create sideset for south boundary
    cubit.cmd(f"group 'boundary_south' add {surfaces_south_str}")
    cubit.cmd("sideset 10 group boundary_south")
    cubit.cmd("sideset 10 name 'boundary_south'")

    # Create sideset for east boundary
    cubit.cmd(f"group 'boundary_east' add {surfaces_east_str}")
    cubit.cmd("sideset 11 group boundary_east")
    cubit.cmd("sideset 11 name 'boundary_east'")

    # Create sideset for north boundary
    cubit.cmd(f"group 'boundary_north' add {surfaces_north_str}")
    cubit.cmd("sideset 12 group boundary_north")
    cubit.cmd("sideset 12 name 'boundary_north'")

    # Create sideset for west boundary
    cubit.cmd(f"group 'boundary_west' add {surfaces_west_str}")
    cubit.cmd("sideset 13 group boundary_west")
    cubit.cmd("sideset 13 name 'boundary_west'")

    # Create sideset for bottom boundary
    cubit.cmd(f"group 'boundary_bottom' add {surfaces_bottom_str}")
    cubit.cmd("sideset 14 group boundary_bottom")
    cubit.cmd("sideset 14 name 'boundary_bottom'")

    # Create sideset for top boundary
    cubit.cmd(f"group 'boundary_top' add {surfaces_top_str}")
    cubit.cmd("sideset 15 group boundary_top")
    cubit.cmd("sideset 15 name 'boundary_top'")


    # Create sideset for fault 'main'
    cubit.cmd(f"group 'fault_main' add {surfaces_fault_main_str}")
    cubit.cmd("sideset 20 group fault_main")
    cubit.cmd("sideset 20 name 'fault_main'")

    # Create sideset for fault 'west'
    cubit.cmd(f"group 'fault_west' add {surfaces_fault_west_str}")
    cubit.cmd("sideset 21 group fault_west")
    cubit.cmd("sideset 21 name 'fault_west'")

    # Create sideset for fault 'east'
    cubit.cmd(f"group 'fault_east' add {surfaces_fault_east_str}")
    cubit.cmd("sideset 22 group fault_east")
    cubit.cmd("sideset 22 name 'fault_east'")


# Write mesh as ExodusII file
cubit.cmd(f"export mesh 'mesh_{CELL}.exo' dimension 3 overwrite")
