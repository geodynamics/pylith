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


# -------------------------------------------------------------------------------------------------
# Mesh definition
# -------------------------------------------------------------------------------------------------
# 2D domain is 20 km x 20 km
# -10 km <= x <= +10 km
# -20 km <= y <= 0
km = 1000.0
LENGTH = 20*km
DEPTH = 20*km

# Magma conduit
CONDUIT_LENGTH = 1.0*km
CONDUIT_DEPTH = 6.0*km

# Magma reservoir
MAJOR_AXIS = 5.0*km
MINOR_AXIS = 1.0*km

# Discretization size in reservoir
DX_RES = 250
DX_BIAS = 1.1

cell = "quad"

# -------------------------------------------------------------------------------------------------
# Start cubit
# -------------------------------------------------------------------------------------------------
setup_cubit()
cubit.reset()


# -------------------------------------------------------------------------------------------------
# Geometry
# -------------------------------------------------------------------------------------------------
# Domain
cubit.cmd(f"create surface rectangle width {LENGTH} height {DEPTH} ZPLANE")
s_domain = get_last_entity("surface")
s_domain.set_entity_name("domain")
v_domain = get_last_entity("volume")
cubit.move(v_domain, [LENGTH/2, -DEPTH/2, 0])

# Chamber
cubit.cmd(f"create surface ellipse major radius {MAJOR_AXIS} minor radius {MINOR_AXIS} on surface domain")
s_chamber = get_last_entity("surface")
s_chamber.set_entity_name("reservoir")
b_chamber = get_last_entity("body")
cubit.move(b_chamber, [LENGTH/2, -DEPTH*(2/3), 0])

# Conduit
cubit.cmd(f"create surface rectangle width {CONDUIT_LENGTH} height {CONDUIT_DEPTH} on surface domain")
b_conduit = get_last_entity("body")
b_conduit.set_entity_name("conduit")
cubit.move(b_conduit, [LENGTH/2, -DEPTH + CONDUIT_DEPTH/2, 0])

# Combine conduit and chamber into reservoir
cubit.unite([b_chamber, b_conduit])

# Imprint to reservoir on domain
cubit.cmd("imprint all")
cubit.cmd("merge all")

# Delete the reservoir
cubit.cmd(f"delete body {b_chamber.id()}")


# -------------------------------------------------------------------------------------------------
# Set names of entities using idless signatures
# -------------------------------------------------------------------------------------------------
# Surfaces
# cubit.get_idless_signature("surface", 4)
cubit.cmd("surface ( at 10000 -16166.7 0 ordinal 1 ) name 's_reservoir'")
# cubit.get_idless_signature("surface", 5)
cubit.cmd("surface ( at 10000 -10000 0 ordinal 1 )name 's_domain'")

# Curves
# cubit.get_idless_signature("curve", 8)
cubit.cmd("curve ( at 10000 -20000 0 ordinal 1 ordered ) name 'c_inflow'")

# cubit.get_idless_signature("curve", 10)
cubit.cmd("curve ( at 10000 -12333.3 0 ordinal 1 ) name 'c_chamber'")

# cubit.get_idless_signature("curve", 11)
cubit.cmd("curve ( at 9500 -17164.2 0 ordinal 1 ) name 'c_conduit_left'")

# cubit.get_idless_signature("curve", 12)
cubit.cmd("curve ( at 10500 -17164.2 0 ordinal 1 ) name 'c_conduit_right'")

# cubit.get_idless_signature("curve", 1)
cubit.cmd("curve ( at 10000 0 0 ordinal 1 ) name 'c_ypos'")

# cubit.get_idless_signature("curve", 17)
cubit.cmd("curve ( at 4750 -20000 0 ordinal 1 ) name 'c_yneg_left'")

# cubit.get_idless_signature("curve", 16)
cubit.cmd("curve ( at 15250 -20000 0 ordinal 1 )  name 'c_yneg_right'")

# cubit.get_idless_signature("curve", 2)
cubit.cmd("curve ( at 0 -10000 0 ordinal 1 ) name 'c_xneg'")

# cubit.get_idless_signature("curve", 4)
cubit.cmd("curve ( at 20000 -10000 0 ordinal 1 ) name 'c_xpos'")


# -------------------------------------------------------------------------------------------------
# Generate the mesh
# -------------------------------------------------------------------------------------------------

# Generate the finite-element mesh.
if cell == "quad":
    cubit.cmd("surface all scheme pave")
else:
    cubit.cmd("surface all scheme trimesh")

# Reset the mesh schemes on curves and the surface sizing function.
cubit.cmd("curve all scheme default")
cubit.cmd("surface all sizing function none")

# Set size on reservoir boundary
cubit.cmd(f"curve c_chamber c_conduit_left c_conduit_right c_inflow size {DX_RES}")

# Use skeleton sizing to increase cell size away from fault.
cubit.cmd(f"surface all sizing function skeleton min_size {DX_RES} max_gradient {DX_BIAS}")

# cubit.cmd("preview mesh surface all")
cubit.cmd("mesh surface all")

# Smooth the mesh to improve mesh quality
cubit.cmd("surface all smooth scheme condition number beta 1.1 cpu 10")
cubit.cmd("smooth surface all")

# -------------------------------------------------------------------------------------------------
# Create blocks for materials and nodesets for boundary conditions.
# -------------------------------------------------------------------------------------------------

# We follow the general approach of creating groups and then creating the nodesets
# from the groups, so that we can apply boolean operations (e.g., union, intersection)
# on the groups to create the desired nodesets.

# Blocks
cubit.cmd("block 1 surface domain")
cubit.cmd("block 1 name 'domain'")

cubit.cmd("block 2 surface reservoir")
cubit.cmd("block 2 name 'reservoir'")

# Nodesets
if True:
    cubit.cmd("group 'boundary_xneg' add node in curve c_xneg")
    cubit.cmd("nodeset 20 group boundary_xneg")
    cubit.cmd("nodeset 20 name 'boundary_xneg'")

    cubit.cmd("group 'boundary_xpos' add node in curve c_xpos")
    cubit.cmd("nodeset 21 group boundary_xpos")
    cubit.cmd("nodeset 21 name 'boundary_xpos'")

    cubit.cmd("group 'boundary_yneg' add node in curve c_yneg_left c_yneg_right c_inflow")
    cubit.cmd("nodeset 22 group boundary_yneg")
    cubit.cmd("nodeset 22 name 'boundary_yneg'")

    cubit.cmd("group 'boundary_ypos' add node in curve c_ypos")
    cubit.cmd("nodeset 23 group boundary_ypos")
    cubit.cmd("nodeset 23 name 'boundary_ypos'")

    cubit.cmd("group 'boundary_flow' add node in curve c_inflow")
    cubit.cmd("nodeset 24 group boundary_flow")
    cubit.cmd("nodeset 24 name 'boundary_flow'")


# Starting in PyLith v5, we will use sidesets instead of nodesets for BCs.
# Sidesets
if False:
    cubit.cmd("group 'boundary_xneg' add curve c_xneg")
    cubit.cmd("sideset 20 group boundary_xneg")
    cubit.cmd("sideset 20 name 'boundary_xneg'")

    cubit.cmd("group 'boundary_xpos' add curve c_xpos")
    cubit.cmd("sideset 21 group boundary_xpos")
    cubit.cmd("sideset 21 name 'boundary_xpos'")

    cubit.cmd("group 'boundary_yneg' add curve c_yneg_left c_yneg_right c_inflow")
    cubit.cmd("sideset 22 group boundary_yneg")
    cubit.cmd("sideset 22 name 'boundary_yneg'")

    cubit.cmd("group 'boundary_ypos' add curve c_ypos")
    cubit.cmd("sideset 23 group boundary_ypos")
    cubit.cmd("sideset 23 name 'boundary_ypos'")

    cubit.cmd("group 'boundary_flow' add curve c_inflow")
    cubit.cmd("sideset 24 group boundary_flow")
    cubit.cmd("sideset 24 name 'boundary_flow'")


# Write mesh as ExodusII file
cubit.cmd(f"export mesh 'mesh_{cell}.exo' dimension 2 overwrite")