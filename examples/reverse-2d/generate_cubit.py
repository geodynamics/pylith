#!python3
"""WARNING: This script only works with Python supplied with Cubit.

To run this script outside the Cubit GUI, you will need to run it like:

PATH_TO_CUBIT/Cubit.app/Contents/Frameworks/Python.framework/Versions/Current/bin/python3 generate_cubit.py --cubit-directory=PATH_TO_CUBIT/Cubit.app/Contents/MacOS

where you replace 'PATH_TO_CUBIT' with the absolute path.
"""

import math

# 2D domain is 200 km x 100 km
# -100 km <= x <= +100 km
# -100 km <= y <= 0
km = 1000.0
BLOCK_WIDTH = 200*km
BLOCK_HEIGHT = 100*km
BLOCK_LENGTH = 10*km

FAULT_WIDTH = 60*km
FAULT_DIP_ANGLE = 30.0 / 180 * math.pi
FAULT_OFFSET = 0.5*FAULT_WIDTH*math.cos(FAULT_DIP_ANGLE)

SPLAY_WIDTH = 15.0*km
SPLAY_DIP_ANGLE = 45.0 / 180.0 * math.pi
SPLAY_OFFSET = 20*km


DX = 2.0*km # Discretization size on fault
BIAS_FACTOR = 1.07 # rate of geometric increase in cell size with distance from the fault
cell = "tri"

# Detect if we are running outside Cubit.
try:
    cubit.reset()
except NameError:
    import argparse
    import pathlib
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument("--cell", action="store", dest="cell", choices=("quad","tri"), default="tri")

    parser.add_argument("--cubit-directory", action="store", dest="cubit_dir", required=True, help="Directory containing cubit executable.")
    args = parser.parse_args()

    # Initialize cubit
    cubit_absdir = pathlib.Path(args.cubit_dir).expanduser().resolve()
    sys.path.append(str(cubit_absdir))
    import cubit
    cubit.init(['cubit','-nojournal'])

    cell = args.cell


cubit.reset()

# -------------------------------------------------------------------------------------------------
# Geometry
# -------------------------------------------------------------------------------------------------

# Create block and then create 2D domain from mid-surface of block.
# Create a brick and move it so fault is centered and upper surface
# is at y=0.
brick = cubit.brick(BLOCK_WIDTH, BLOCK_HEIGHT, BLOCK_LENGTH)
cubit.move(brick, (-FAULT_OFFSET, -0.5*BLOCK_HEIGHT, 0))

# Create a midsurface from front and back surfaces.
# surface 1 name "surf_front"
# surface 2 name "surf_back"
cubit.cmd("surface  ( at -25980.8 -50000 5000 ordinal 1 ordered )  name 'surf_front'")
cubit.cmd("surface  ( at -25980.8 -50000 -5000 ordinal 1 ordered )  name 'surf_back'")
cubit.cmd(f"create midsurface volume {brick.id()} surface surf_front surf_back")
surf_id = cubit.get_last_id("surface")

# Delete the initial volume now we have the midsurface.
cubit.cmd(f"delete volume {brick.id()}")


# Create fault and splay surfaces

# Create fault (yz plane) at x = 0.0
cubit.cmd(f"split surface {surf_id} across location position 0 0 0 location position {-BLOCK_HEIGHT/math.tan(FAULT_DIP_ANGLE)} {-BLOCK_HEIGHT} 0")

# split curve 17 at position {-faultWidth*cosd(fault1DipAngle)} {-faultWidth*sind(faultDipAngle)} 0
cubit.cmd(f"split curve  ( at -62990.4 -36367.5 0 ordinal 1 ordered )  at position {-FAULT_WIDTH*math.cos(FAULT_DIP_ANGLE)} {-FAULT_WIDTH*math.sin(FAULT_DIP_ANGLE)} 0")

# split surface 8 across location position -{splayOffset} 0 0 location position {-BLOCK_HEIGHT/tand(splayDipAngle)} {-BLOCK_HEIGHT} 0
cubit.cmd(f"split surface  ( at -62990.4 -36367.5 0 ordinal 1 ordered )  across location position {-SPLAY_OFFSET} 0 0 location position {-BLOCK_HEIGHT/math.tan(SPLAY_DIP_ANGLE)} {-BLOCK_HEIGHT} 0")


# Name surfaces

# Material on the footwall.
# surface 9 name "surface_slab"
cubit.cmd("surface  ( at -25980.8 -50000 0 ordinal 1 ordered )  name 'surface_slab'")

# Material on the hanging wall.
# surface 10 name 'surface_crust'
cubit.cmd("surface  ( at -72990.4 -36367.5 0 ordinal 1 ordered )  name 'surface_crust'")

# Material between faul and splay.
# surface 11 name 'surface_wedge'
cubit.cmd("surface  ( at -18583.2 -10729 0 ordinal 1 ordered )  name 'surface_wedge'")


# Name curves

# curve 26 name 'c_fault_lower'
cubit.cmd("curve  ( at -44564 -25729 0 ordinal 1 ordered )  name 'c_fault_lower'")

# curve 27 name 'c_fault_upper'
cubit.cmd("curve  ( at -18583.2 -10729 0 ordinal 1 ordered )  name 'c_fault_upper'")

# curve 24 name 'c_splay'
cubit.cmd("curve  ( at -28583.2 -10729 0 ordinal 1 ordered )  name 'c_splay'")

# curve 21 name 'c_ypos_fw'
cubit.cmd("curve  ( at 37009.6 0 0 ordinal 1 ordered )  name 'c_ypos_fw'")

# curve 28 name 'c_ypos_w'
cubit.cmd("curve  ( at -10000 0 0 ordinal 1 ordered )  name 'c_ypos_w'")

# curve 25 name 'c_ypos_hw'
cubit.cmd("curve  ( at -72990.4 0 0 ordinal 1 ordered )  name 'c_ypos_hw'")

# curve 14 name 'c_xpos'
cubit.cmd("curve  ( at 74019.2 -50000 0 ordinal 1 ordered )  name 'c_xpos'")

# curve 20 name 'c_xneg_fw'
cubit.cmd("curve  ( at -125981 -86367.5 0 ordinal 1 ordered )  name 'c_xneg_fw'")

# curve 19 name 'c_xneg_hw'
cubit.cmd("curve  ( at -125981 -36367.5 0 ordinal 1 ordered )  name 'c_xneg_hw'")

# curve 13 name 'c_yneg'
cubit.cmd("curve  ( at -25980.8 -100000 0 ordinal 1 ordered )  name 'c_yneg'")

# curve 23 name 'c_fault_ext'
cubit.cmd("curve  ( at -88971.1 -51367.5 0 ordinal 1 ordered )  name 'c_fault_ext'")


# Name vertices

# vertex 15 name 'v_fault_bot'
cubit.cmd("vertex  ( at -51961.5 -30000 0 ordinal 1 ordered )  name 'v_fault_bot'")

# vertex 13 name 'v_fault_top'
cubit.cmd("vertex  ( at 0 0 0 ordinal 1 ordered )  name 'v_fault_top'")

# vertex 17 name 'v_splay_bot'
cubit.cmd("vertex  ( at -37166.4 -21458.1 0 ordinal 1 ordered )  name 'v_splay_bot'")

# vertex 16 name 'v_splay_top'
cubit.cmd("vertex  ( at -20000 -3.63798e-12 0 ordinal 1 ordered )  name 'v_splay_top'")

# vertex 14 name 'v_fault_xneg'
cubit.cmd("vertex  ( at -125981 -72735 0 ordinal 1 ordered )  name 'v_fault_xneg'")

# vertex 11 name 'v_ypos_xpos'
cubit.cmd("vertex  ( at 74019.2 0 0 ordinal 1 ordered )  name 'v_ypos_xpos'")

# vertex 12 name 'v_ypos_xneg'
cubit.cmd("vertex  ( at -125981 0 0 ordinal 1 ordered )  name 'v_ypos_xneg'")

# vertex 9 name 'v_yneg_xpos'
cubit.cmd("vertex  ( at 74019.2 -100000 0 ordinal 1 ordered )  name 'v_yneg_xpos'")

# vertex 10 name 'v_yneg_xneg'
cubit.cmd("vertex  ( at -125981 -100000 0 ordinal 1 ordered )  name 'v_yneg_xneg'")


# -------------------------------------------------------------------------------------------------
# Generate the mesh
# -------------------------------------------------------------------------------------------------

# Generate the finite-element mesh.
if cell == "quad":
    cubit.cmd("surface all scheme pave")
else:
    cubit.cmd("surface all scheme trimesh")

# Set sizes to create mesh with cell size than increases at a geometric rate with distance from the fault.
def compute_dx_curve_end(dx_start, curve_length):
    """Compute cell size at end of curve given cell size at start of curve."""
    return dx_start * BIAS_FACTOR**math.ceil(math.log(1-curve_length/dx_start*(1-BIAS_FACTOR))/math.log(BIAS_FACTOR))

# Compute sizes at curve endpoints
# ----------------------------------------------------------------------

# dxA - size at v_top_xpos
curve_id = cubit.get_id_from_name("c_ypos_fw")
dx_A = compute_dx_curve_end(dx_start=DX, curve_length=cubit.get_curve_length(curve_id))

# dxB - size at v_top_xneg
curve_id = cubit.get_id_from_name("c_ypos_hw")
dx_B = compute_dx_curve_end(dx_start=DX, curve_length=cubit.get_curve_length(curve_id))

# dxC - size at v_fault_xneg
curve_id = cubit.get_id_from_name("c_fault_ext")
dx_C = compute_dx_curve_end(dx_start=DX, curve_length=cubit.get_curve_length(curve_id))

# dxD - size at v_bot_xpos
curve_id = cubit.get_id_from_name("c_xpos")
dx_D = compute_dx_curve_end(dx_start=dx_A, curve_length=cubit.get_curve_length(curve_id))


# Reset sizes
cubit.cmd("curve all scheme default")
cubit.cmd("surface all sizing function none")

# Set size on faults
cubit.cmd(f"curve c_fault_upper size {DX}")
cubit.cmd(f"curve c_fault_lower size {DX}")
cubit.cmd(f"curve c_splay size {DX}")
cubit.cmd(f"curve c_ypos_w size {DX}")

# Set bias on curves extending from faults
cubit.cmd(f"curve c_ypos_fw scheme bias fine size {DX} factor {BIAS_FACTOR} start vertex v_fault_top")
cubit.cmd(f"curve c_ypos_hw scheme bias fine size {DX} factor {BIAS_FACTOR} start vertex v_splay_top")
cubit.cmd(f"curve c_fault_ext scheme bias fine size {DX} factor {BIAS_FACTOR} start vertex v_fault_bot")

cubit.cmd(f"curve c_yneg size {dx_D}")

# A to D
cubit.cmd(f"curve c_xpos scheme bias fine size {dx_A} coarse size {dx_D} start vertex v_ypos_xpos")

# C to B
cubit.cmd(f"curve c_xneg_hw scheme bias fine size {dx_C} coarse size {dx_B} start vertex v_fault_xneg")

# B to D
cubit.cmd(f"curve c_xneg_fw scheme bias fine size {dx_B} coarse size {dx_D} start vertex v_fault_xneg")

cubit.cmd(f"surface all sizing function type bias start curve c_fault_upper c_fault_lower c_splay c_ypos_w factor {BIAS_FACTOR}")


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
cubit.cmd("block 1 surface surface_slab")
cubit.cmd("block 1 name 'surface_slab'")

cubit.cmd("block 2 surface surface_crust")
cubit.cmd("block 2 name 'surface_crust'")

cubit.cmd("block 3 surface surface_wedge")
cubit.cmd("block 3 name 'surface_wedge'")


# Nodesets

# Create nodeset for fault
cubit.cmd("group 'fault' add node in c_fault_upper c_fault_lower")
cubit.cmd("nodeset 10 group fault")
cubit.cmd("nodeset 10 name 'fault'")

# Create nodeset for fault edge
cubit.cmd("group 'fault_end' add node in vertex v_fault_bot")
cubit.cmd("nodeset 11 group fault_end")
cubit.cmd("nodeset 11 name 'fault_end'")

# Create nodeset for splay
cubit.cmd("group 'splay' add node in c_splay")
cubit.cmd("nodeset 12 group splay")
cubit.cmd("nodeset 12 name 'splay'")

# Create nodeset for splay edge
cubit.cmd("group 'splay_end' add node in vertex v_splay_bot")
cubit.cmd("nodeset 13 group splay_end")
cubit.cmd("nodeset 13 name 'splay_end'")

# Create nodeset for +x edge
cubit.cmd("group 'boundary_xpos' add node in curve c_xpos")
cubit.cmd("nodeset 20 group boundary_xpos")
cubit.cmd("nodeset 20 name 'boundary_xpos'")

# Create nodeset for -x edge
cubit.cmd("group 'boundary_xneg' add node in curve c_xneg_hw")
cubit.cmd("group 'boundary_xneg' add node in curve c_xneg_fw")
cubit.cmd("nodeset 21 group boundary_xneg")
cubit.cmd("nodeset 21 name 'boundary_xneg'")

# Create nodeset for +y edge
cubit.cmd("group 'boundary_ypos' add node in curve c_ypos_fw c_ypos_hw c_ypos_w")
cubit.cmd("nodeset 22 group boundary_ypos")
cubit.cmd("nodeset 22 name 'boundary_ypos'")

# Create nodeset for -y edge
cubit.cmd("group 'boundary_yneg' add node in curve c_yneg")
cubit.cmd("nodeset 23 group boundary_yneg")
cubit.cmd("nodeset 23 name 'boundary_yneg'")


# Starting in PyLith v5, we will use sidesets instead of nodesets for BCs.
# Sidesets
#
# Create sideset for fault
#cubit.cmd("group 'fault' add c_fault_upper c_fault_lower")
#cubit.cmd("sideset 10 group fault")
#cubit.cmd("sideset 10 name 'fault'")
#
# Create sideset for splay
#cubit.cmd("group 'splay' add c_splay")
#cubit.cmd("sideset 12 group splay")
#cubit.cmd("sideset 12 name 'splay'")
#
# Create sideset for +x edge
#cubit.cmd("group 'boundary_xpos' add curve c_xpos")
#cubit.cmd("sideset 20 group boundary_xpos")
#cubit.cmd("sideset 20 name 'boundary_xpos'")
#
# Create sideset for -x edge
#cubit.cmd("group 'boundary_xneg' add curve c_xneg_hw")
#cubit.cmd("group 'boundary_xneg' add curve c_xneg_fw")
#cubit.cmd("sideset 21 group boundary_xneg")
#cubit.cmd("sideset 21 name 'boundary_xneg'")
#
# Create sideset for +y edge
#cubit.cmd("group 'boundary_ypos' add curve c_ypos_fw c_ypos_hw c_ypos_w")
#cubit.cmd("sideset 22 group boundary_ypos")
#cubit.cmd("sideset 22 name 'boundary_ypos'")
#
# Create sideset for -y edge
#cubit.cmd("group 'boundary_yneg' add curve c_yneg")
#cubit.cmd("sideset 23 group boundary_yneg")
#cubit.cmd("sideset 23 name 'boundary_yneg'")


# Write mesh as ExodusII file
cubit.cmd(f"export mesh 'mesh_{cell}.exo' dimension 2 overwrite")
