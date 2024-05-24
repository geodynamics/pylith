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


def entity_from_name(name: str, entity_type: str):
    get_type = {
        "vertex": cubit.vertex,
        "curve": cubit.curve,
        "surface": cubit.surface,
        "volume": cubit.volume,
        "body": cubit.body,
    }
    return get_type[entity_type](cubit.get_id_from_name(name))


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
# We build the domain by forming curves from points and the surface from curves. We
# imprint the fault traces on the surface. With this workflow we name objects as we
# create them.

import numpy

# Locations defining domain
DOMAIN_W = 410000.00
DOMAIN_E = 490000.00
DOMAIN_S = 3910000.00
DOMAIN_N = 3990000.00

# Files with coordinates of fault traces in UTM zone 11
FILENAME_MAINTRACE = "faulttrace_main_utm.txt"
FILENAME_WESTTRACE = "faulttrace_west_utm.txt"
FILENAME_STRANDETRACE = "faulttrace_east_utm.txt"

# Discretization size on faults and bias (rate of cell size increase away from fault)
DX_FAULT = 1.0e+3
DX_BIAS = 1.07

CELL = "tri"

# -------------------------------------------------------------------------------------------------
# Start cubit
# -------------------------------------------------------------------------------------------------
setup_cubit()
cubit.reset()


# -------------------------------------------------------------------------------------------------
# Geometry
# -------------------------------------------------------------------------------------------------
# Create points for domain
pSW = cubit.create_vertex(DOMAIN_W, DOMAIN_S, 0)
pSE = cubit.create_vertex(DOMAIN_E, DOMAIN_S, 0)
pNE = cubit.create_vertex(DOMAIN_E, DOMAIN_N, 0)
pNW = cubit.create_vertex(DOMAIN_W, DOMAIN_N, 0)

# Create curves for domain
c_south = cubit.create_curve(pSW, pSE); c_south.set_entity_name("c_south")
c_east = cubit.create_curve(pSE, pNE); c_east.set_entity_name("c_east")
c_north = cubit.create_curve(pNE, pNW); c_north.set_entity_name("c_north")
c_west = cubit.create_curve(pNW, pSW); c_west.set_entity_name("c_west")

# Create surface for domain
s_domain = cubit.create_surface([c_south, c_east, c_north, c_west]); s_domain.set_entity_name("s_domain")


# Create points for main fault
points_fault_main = create_points_from_file(FILENAME_MAINTRACE)
cubit.cmd(f"create curve spline vertex {points_fault_main[0].id()} to {points_fault_main[-1].id()}")
c_fault_main = get_last_entity("curve")
cubit.cmd(f"curve {c_fault_main.id()} name 'c_fault_main'")

# Create points for west fault strand
points_fault_west = create_points_from_file(FILENAME_WESTTRACE)
cubit.cmd(f"create curve spline vertex {points_fault_west[0].id()} to {points_fault_west[-1].id()}")
c_fault_west = get_last_entity("curve")
cubit.cmd(f"curve {c_fault_west.id()} name 'c_fault_west'")

# Create points for east fault strand
points_fault_east = create_points_from_file(FILENAME_STRANDETRACE)
cubit.cmd(f"create curve spline vertex {points_fault_east[0].id()} to {points_fault_east[-1].id()}")
c_fault_east = get_last_entity("curve")
cubit.cmd(f"curve {c_fault_east.id()} name 'c_fault_east'")

# Imprint/merge
cubit.cmd("delete vertex all")
cubit.cmd(f"imprint surface {s_domain.id()} with curve {c_fault_main.id()}")
cubit.cmd(f"imprint surface {s_domain.id()} with curve {c_fault_west.id()}")
cubit.cmd(f"imprint surface {s_domain.id()} with curve {c_fault_east.id()}")
cubit.cmd("merge all")

# Merging destroys id values but not names, so use names (when possible) to get update entities
fault_west_ends = entity_from_name("c_fault_west", "curve").vertices()
fault_east_ends = entity_from_name("c_fault_east", "curve").vertices()

# Manually find ids using GUI for main fault split by imprint, then replace id with idless signature
#idless = cubit.get_idless_signature("curve", 9)
idless = "Curve ( at 447417 3.95436e+06 0 ordinal 1 )"
cubit.cmd(f"curve {idless} name 'c_fault_main_north'")
#idless = cubit.get_idless_signature("curve", 10)
idless = "Curve ( at 460642 3.94255e+06 0 ordinal 1 )"
cubit.cmd(f"curve {idless} name 'c_fault_main_south'")
fault_main_ends = [cubit.vertex(9), cubit.vertex(15)]


# -------------------------------------------------------------------------------------------------
# Generate the mesh
# -------------------------------------------------------------------------------------------------
cubit.cmd("surface all scheme trimesh")

# Reset sizes
cubit.cmd("curve all scheme default")
cubit.cmd("surface all sizing function none")

# Set size on faults
cubit.cmd(f"curve c_fault_main_north c_fault_main_south size {DX_FAULT}")
cubit.cmd(f"curve c_fault_west size {DX_FAULT}")
cubit.cmd(f"curve c_fault_west size {DX_FAULT}")

cubit.cmd(f"surface all sizing function skeleton min_size {DX_FAULT} max_gradient {DX_BIAS}")

# cubit.cmd("preview mesh surface all")
cubit.cmd("mesh surface all")

# Smooth the mesh to improve mesh quality
cubit.cmd("surface all smooth scheme condition number beta 1.3 cpu 10")
cubit.cmd("smooth surface all")


# -------------------------------------------------------------------------------------------------
# Create blocks for materials and nodesets for boundary conditions.
# -------------------------------------------------------------------------------------------------

# We follow the general approach of creating groups and then creating the nodesets
# from the groups, so that we can apply boolean operations (e.g., union, intersection)
# on the groups to create the desired nodesets.

# Blocks
cubit.cmd(f"block 1 surface {s_domain.id()}")
cubit.cmd("block 1 name 'domain'")


# Nodesets

if True:
    # Create nodeset for south boundary
    cubit.cmd(f"group 'boundary_south' add node in curve {c_south.id()}")
    cubit.cmd("nodeset 10 group boundary_south")
    cubit.cmd("nodeset 10 name 'boundary_south'")

    # Create nodeset for east boundary
    cubit.cmd(f"group 'boundary_east' add node in curve {c_east.id()}")
    cubit.cmd("nodeset 11 group boundary_east")
    cubit.cmd("nodeset 11 name 'boundary_east'")

    # Create nodeset for north boundary
    cubit.cmd(f"group 'boundary_north' add node in curve {c_north.id()}")
    cubit.cmd("nodeset 12 group boundary_north")
    cubit.cmd("nodeset 12 name 'boundary_north'")

    # Create nodeset for west boundary
    cubit.cmd(f"group 'boundary_west' add node in curve {c_west.id()}")
    cubit.cmd("nodeset 13 group boundary_west")
    cubit.cmd("nodeset 13 name 'boundary_west'")


    # Create nodeset for fault 'main'
    cubit.cmd("group 'fault_main' add node in curve c_fault_main_north c_fault_main_south")
    cubit.cmd("nodeset 20 group fault_main")
    cubit.cmd("nodeset 20 name 'fault_main'")

    # Create nodeset for fault 'main' ends
    cubit.cmd(f"group 'fault_main_ends' add node in vertex {ids_from_entities(fault_main_ends)}")
    cubit.cmd("nodeset 30 group fault_main_ends")
    cubit.cmd("nodeset 30 name 'fault_main_ends'")

    # Create nodeset for fault 'west'
    cubit.cmd("group 'fault_west' add node in curve c_fault_west")
    cubit.cmd("nodeset 21 group fault_west")
    cubit.cmd("nodeset 21 name 'fault_west'")

    # Create nodeset for fault 'west' ends
    cubit.cmd(f"group 'fault_west_ends' add node in vertex {ids_from_entities(fault_west_ends)}")
    cubit.cmd("nodeset 31 group fault_west_ends")
    cubit.cmd("nodeset 31 name 'fault_west_ends'")

    # Create nodeset for fault 'east'
    cubit.cmd("group 'fault_east' add node in curve c_fault_east")
    cubit.cmd("nodeset 22 group fault_east")
    cubit.cmd("nodeset 22 name 'fault_east'")

    # Create nodeset for fault 'east' ends
    cubit.cmd(f"group 'fault_east_ends' add node in vertex {ids_from_entities(fault_east_ends)}")
    cubit.cmd("nodeset 32 group fault_east_ends")
    cubit.cmd("nodeset 32 name 'fault_east_ends'")


# Starting in PyLith v5, we will use sidesets instead of nodesets for BCs.
# Sidesets

if False:
    # Create sideset for south boundary
    cubit.cmd(f"group 'boundary_south' add curve {c_south.id()}")
    cubit.cmd("sideset 10 group boundary_south")
    cubit.cmd("sideset 10 name 'boundary_south'")

    # Create sideset for east boundary
    cubit.cmd(f"group 'boundary_east' add curve {c_east.id()}")
    cubit.cmd("sideset 11 group boundary_east")
    cubit.cmd("sideset 11 name 'boundary_east'")

    # Create sideset for north boundary
    cubit.cmd(f"group 'boundary_north' add curve {c_north.id()}")
    cubit.cmd("sideset 12 group boundary_north")
    cubit.cmd("sideset 12 name 'boundary_north'")

    # Create sideset for west boundary
    cubit.cmd(f"group 'boundary_west' add curve {c_west.id()}")
    cubit.cmd("sideset 13 group boundary_west")
    cubit.cmd("sideset 13 name 'boundary_west'")


    # Create sideset for fault 'main'
    cubit.cmd("group 'fault_main' add curve c_fault_main_north c_fault_main_south")
    cubit.cmd("sideset 20 group fault_main")
    cubit.cmd("sideset 20 name 'fault_main'")

    # Create sideset for fault 'west'
    cubit.cmd("group 'fault_west' add curve c_fault_west")
    cubit.cmd("sideset 21 group fault_west")
    cubit.cmd("sideset 21 name 'fault_west'")

    # Create sideset for fault 'east'
    cubit.cmd("group 'fault_east' add curve c_fault_east")
    cubit.cmd("sideset 22 group fault_east")
    cubit.cmd("sideset 22 name 'fault_east'")


# Write mesh as ExodusII file
cubit.cmd(f"export mesh 'mesh_{CELL}.exo' dimension 2 overwrite")
