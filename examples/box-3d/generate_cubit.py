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
        parser.add_argument("--cell", action="store", dest="cell", choices=("hex","tet"), default="hex")

        parser = argparse.ArgumentParser()
        parser.add_argument("--cubit-directory", action="store", dest="cubit_dir", required=True, help="Directory containing cubit executable.")
        args = parser.parse_args()

        # Initialize cubit
        cubit_absdir = pathlib.Path(args.cubit_dir).expanduser().resolve()
        sys.path.append(str(cubit_absdir))
        import cubit
        cubit.init(['cubit','-nojournal'])

# -------------------------------------------------------------------------------------------------
# Mesh definition
# -------------------------------------------------------------------------------------------------
km = 1000.0
DOMAIN_X = DOMAIN_Y = 12.0*km
DOMAIN_Z = 9.0*km
DX = 3.0*km
cell = "hex"

# -------------------------------------------------------------------------------------------------
# Start cubit
# -------------------------------------------------------------------------------------------------
setup_cubit()
cubit.reset()

# -------------------------------------------------------------------------------------------------
# Geometry
# -------------------------------------------------------------------------------------------------
# Create block.
#
# -6 km <= x <= 6 km
# -6 km <= y <= 6 km
# -9 km <= z <= 0 km
brick = cubit.brick(DOMAIN_X, DOMAIN_Y, DOMAIN_Z)
cubit.move(brick, (0, 0, -0.5*DOMAIN_Z))

# Create a mapping of surface 'id' to name, so we can give the surfaces meaningful names.
# The surface 'id' is determined by Cubit, so we examined the 'id' of each surface to determine
# which 'id' corresponds to each name.
surfaces = { \
    1: "surface_zpos", \
    2: "surface_zneg", \
    3: "surface_yneg", \
    4: "surface_xneg", \
    5: "surface_ypos", \
    6: "surface_xpos", \
}
for s_id, surface_name in surfaces.items():
    cubit.cmd(f"surface {s_id} name '{surface_name}'")

# Generate the finite-element mesh.
cubit.cmd(f"volume all size {DX}")
if cell == "hex":
    cubit.cmd("volume all scheme map")
else:
    cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# Create blocks for materials and nodesets for boundary conditions.

# We follow the general approach of creating groups and then creating the nodesets
# from the groups, so that we can apply boolean operations (e.g., union, intersection)
# on the groups to create the desired nodesets.

# Blocks
cubit.cmd(f"block 1 volume {brick.id()}")
cubit.cmd("block 1 name 'elastic'")

# Nodesets
#nodeset_id = 20
#for surface_name in surfaces.values():
#    group_name = surface_name.replace("surface", "boundary")
#    nodeset_name = group_name
#    cubit.cmd(f"group '{group_name}' add node in surface {surface_name}")
#    cubit.cmd(f"nodeset {nodeset_id} group {group_name}")
#    cubit.cmd(f"nodeset {nodeset_id} name '{nodeset_name}'")
#    nodeset_id += 1


# Starting in PyLith v5, we will use sidesets instead of nodesets for BCs.
# Sidesets
sideset_id = 20
for surface_name in surfaces.values():
    group_name = surface_name.replace("surface", "boundary")
    sideset_name = group_name
    cubit.cmd(f"group '{group_name}' add surface {surface_name}")
    cubit.cmd(f"sideset {sideset_id} group {group_name}")
    cubit.cmd(f"sideset {sideset_id} name '{sideset_name}'")
    sideset_id += 1

# Write mesh as ExodusII file
cubit.cmd(f"export mesh 'mesh_{cell}.exo' dimension 3 overwrite")
