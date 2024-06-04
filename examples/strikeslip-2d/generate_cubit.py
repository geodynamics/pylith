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


# -------------------------------------------------------------------------------------------------
# Mesh definition
# -------------------------------------------------------------------------------------------------
# 2D domain is 100 km x 150 km
# -50 km <= x <= +50 km
# -75 km <= y <= +75 km
km = 1000.0
BLOCK_WIDTH = 100*km
BLOCK_LENGTH = 150*km
BLOCK_HEIGHT = 10*km

DX_FAULT = 8.0*km # Discretization size on fault
DX_BIAS = 1.1 # rate of geometric increase in cell size with distance from the fault
cell = "quad"

# -------------------------------------------------------------------------------------------------
# Start cubit
# -------------------------------------------------------------------------------------------------
setup_cubit()
cubit.reset()

import math


# -------------------------------------------------------------------------------------------------
# Geometry
# -------------------------------------------------------------------------------------------------
# Create block and then create 2D domain from mid-surface of block.
brick = cubit.brick(BLOCK_WIDTH, BLOCK_LENGTH, BLOCK_HEIGHT)

cubit.cmd(f"surface ( at 0 0 {0.5*BLOCK_HEIGHT} ordinal 1 ordered ) name 'surf_front'")
cubit.cmd(f"surface ( at 0 0 {-0.5*BLOCK_HEIGHT} ordinal 1 ordered ) name 'surf_back'")
cubit.cmd(f"create midsurface volume {brick.id()} surface surf_front surf_back")
domain_id = cubit.get_last_id("surface")
cubit.cmd(f"delete volume {brick.id()}")

# Create fault (yz plane) at x = 0.0
cubit.cmd(f"split surface {domain_id} across location position 0 {-0.5*BLOCK_LENGTH} 0 location position 0 {+0.5*BLOCK_LENGTH} 0")
cubit.cmd("curve ( at 0 0 0 ordinal 1 ordered ) name 'fault'")

# Name surfaces split by fault interface
cubit.cmd(f"surface  ( at {+0.25*BLOCK_WIDTH} 0 0 ordinal 1 ordered )  name 'surface_xpos'")
cubit.cmd(f"surface  ( at {-0.25*BLOCK_WIDTH} 0 0 ordinal 1 ordered )  name 'surface_xneg'")

# Name curves
cubit.cmd(f"curve ( at {+0.25*BLOCK_WIDTH} {+0.5*BLOCK_LENGTH} 0 ordinal 1 ordered ) name 'c_ypos_xpos'")
cubit.cmd(f"curve ( at {-0.25*BLOCK_WIDTH} {+0.5*BLOCK_LENGTH} 0 ordinal 1 ordered ) name 'c_ypos_xneg'")
cubit.cmd(f"curve ( at {+0.25*BLOCK_WIDTH} {-0.5*BLOCK_LENGTH} 0 ordinal 1 ordered ) name 'c_yneg_xpos'")
cubit.cmd(f"curve ( at {-0.25*BLOCK_WIDTH} {-0.5*BLOCK_LENGTH} 0 ordinal 1 ordered ) name 'c_yneg_xneg'")
cubit.cmd(f"curve ( at {+0.5*BLOCK_WIDTH} 0 0 ordinal 1 ordered ) name 'c_xpos'")
cubit.cmd(f"curve ( at {-0.5*BLOCK_WIDTH} 0 0 ordinal 1 ordered ) name 'c_xneg'")

# Name vertices
cubit.cmd(f"vertex ( at 0 {+0.5*BLOCK_LENGTH} 0 ordinal 1 ordered ) name 'v_fault_ypos'")
cubit.cmd(f"vertex ( at 0 {-0.5*BLOCK_LENGTH} 0 ordinal 1 ordered ) name 'v_fault_yneg'")
cubit.cmd(f"vertex ( at {+0.5*BLOCK_WIDTH} {+0.5*BLOCK_LENGTH} 0 ordinal 1 ordered ) name 'v_ypos_xpos'")
cubit.cmd(f"vertex ( at {-0.5*BLOCK_WIDTH} {+0.5*BLOCK_LENGTH} 0 ordinal 1 ordered ) name 'v_ypos_xneg'")
cubit.cmd(f"vertex ( at {+0.5*BLOCK_WIDTH} {-0.5*BLOCK_LENGTH} 0 ordinal 1 ordered ) name 'v_yneg_xpos'")
cubit.cmd(f"vertex ( at {-0.5*BLOCK_WIDTH} {-0.5*BLOCK_LENGTH} 0 ordinal 1 ordered ) name 'v_yneg_xneg'")

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

# Set size on faults
cubit.cmd(f"curve fault size {DX_FAULT}")

# Use skeleton sizing to increase cell size away from fault.
cubit.cmd(f"surface all sizing function skeleton min_size {DX_FAULT} max_gradient {DX_BIAS}")

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
cubit.cmd("block 1 surface surface_xneg")
cubit.cmd("block 1 name 'elastic_xneg'")

cubit.cmd("block 2 surface surface_xpos")
cubit.cmd("block 2 name 'elastic_xpos'")


# Nodesets
if True:
    cubit.cmd("group 'fault_g' add curve fault")
    cubit.cmd("nodeset 10 group fault_g")
    cubit.cmd("nodeset 10 name 'fault'")

    cubit.cmd("group 'boundary_xpos' add node in curve c_xpos")
    cubit.cmd("nodeset 21 group boundary_xpos")
    cubit.cmd("nodeset 21 name 'boundary_xpos'")

    cubit.cmd("group 'boundary_xneg' add node in curve c_xneg")
    cubit.cmd("nodeset 22 group boundary_xneg")
    cubit.cmd("nodeset 22 name 'boundary_xneg'")

    cubit.cmd("group 'boundary_ypos' add node in curve c_ypos_xneg")
    cubit.cmd("group 'boundary_ypos' add node in curve c_ypos_xpos")
    cubit.cmd("nodeset 23 group boundary_ypos")
    cubit.cmd("nodeset 23 name 'boundary_ypos'")

    cubit.cmd("group 'boundary_yneg' add node in curve c_yneg_xneg")
    cubit.cmd("group 'boundary_yneg' add node in curve c_yneg_xpos")
    cubit.cmd("nodeset 24 group boundary_yneg")
    cubit.cmd("nodeset 24 name 'boundary_yneg'")


# Starting in PyLith v5, we will use sidesets instead of nodesets for BCs.
# Sidesets
if False:
    cubit.cmd("group 'fault_g' add curve fault")
    cubit.cmd("sideset 10 group fault_g")
    cubit.cmd("sideset 10 name 'fault'")

    cubit.cmd("group 'boundary_xpos' add curve c_xpos")
    cubit.cmd("sideset 21 group boundary_xpos")
    cubit.cmd("sideset 21 name 'boundary_xpos'")

    cubit.cmd("group 'boundary_xneg' add curve c_xneg")
    cubit.cmd("sideset 22 group boundary_xneg")
    cubit.cmd("sideset 22 name 'boundary_xneg'")

    cubit.cmd("group 'boundary_ypos' add curve c_ypos_xneg")
    cubit.cmd("group 'boundary_ypos' add curve c_ypos_xpos")
    cubit.cmd("sideset 23 group boundary_ypos")
    cubit.cmd("sideset 23 name 'boundary_ypos'")

    cubit.cmd("group 'boundary_yneg' add curve c_yneg_xneg")
    cubit.cmd("group 'boundary_yneg' add curve c_yneg_xpos")
    cubit.cmd("sideset 24 group boundary_yneg")
    cubit.cmd("sideset 24 name 'boundary_yneg'")


# Write mesh as ExodusII file
cubit.cmd(f"export mesh 'mesh_{cell}.exo' dimension 2 overwrite")
