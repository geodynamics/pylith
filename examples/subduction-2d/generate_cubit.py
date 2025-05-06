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
# 2D subduction zone example based on the 2011 M9.0 Tohoku earthquake.
km = 1000.0
DX_FAULT = 5.0*km # Discretization size on fault
DX_BIAS = 1.07 # rate of geometric increase in cell size with distance from the fault
CELL = "tri"


# -------------------------------------------------------------------------------------------------
# Start cubit
# -------------------------------------------------------------------------------------------------
setup_cubit()
cubit.reset()

import math


# -------------------------------------------------------------------------------------------------
# Geometry
# -------------------------------------------------------------------------------------------------
# The steps in constructing the geometry are:
# (1) Create points
# (2) Connect the points into spline curves
# (3) Split the splines at intersections to form bounding curves
# (4) Form surfaces from the bounding curves
#
# Points have been projected from longitude/latitude into a local
# transverse Mercator projection. PyLith uses the Proj.4 library
# for geographic projections. The proj parameters are:
#
#   +proj=tmerc +datum=WGS84 +lon_0=142.0 +lat_0=38.0 +k=0.9996
#
# so that the local origin is at a longitude of 142.0 degrees (WGS84)
# and a latitude of 38.0 degrees (WGS84).

# Create block and then create 2D domain from mid-surface of block.
# Create a brick and move it so fault is centered and upper surface
# is at y=0.

# Topography/bathymetry (points from Google Earth)
#
# Save ids of important points in APREPRO variables as they are created
# for use in other commands. We save the entity id in a variable rather
# than name the vertex because this allows us to still use "idA to idB"
# to select multiple vertices.
#
p_topo_west = cubit.create_vertex(-600.0*km, -2.0*km, 0.0)
cubit.create_vertex(-439.1*km, -0.3*km, 0.0)
cubit.create_vertex(-351.2*km, -0.8*km, 0.0)
cubit.create_vertex(-263.4*km,  0.0*km, 0.0)
cubit.create_vertex(-175.6*km,  0.4*km, 0.0)
cubit.create_vertex( -87.7*km,  0.0*km, 0.0)
cubit.create_vertex(   0.0*km, -0.4*km, 0.0)
cubit.create_vertex(  87.7*km, -3.0*km, 0.0)
p_topo_trench = cubit.create_vertex(165.6*km, -6.0*km, 0.0)
cubit.create_vertex( 263.4*km, -5.40*km, 0.0)
cubit.create_vertex( 351.2*km, -5.40*km, 0.0)
cubit.create_vertex( 439.1*km, -5.40*km, 0.0)
p_topo_east = cubit.create_vertex( 600.0*km, -5.70*km, 0.0)

cubit.cmd(f"create curve spline vertex {p_topo_west.id()} to {p_topo_east.id()}")
cubit.cmd(f"curve {cubit.get_last_id('curve')} name 'c_topo'")

# Top of slab
# Hayes and Wald, 2009
# http://earthquake.usgs.gov/research/data/slab
p_slabtop_west = cubit.create_vertex(-422.4*km, -240.00*km, 0.0)
cubit.create_vertex(-331.0*km, -180.00*km, 0.0)
cubit.create_vertex(-261.6*km, -140.00*km, 0.0)
cubit.create_vertex(-223.9*km, -120.00*km, 0.0)
cubit.create_vertex(-182.6*km, -100.00*km, 0.0)
cubit.create_vertex(-134.3*km, -80.00*km, 0.0)
p_slabtop_coseismic = cubit.create_vertex( -74.6*km, -60.00*km, 0.0)
p_slabtop_moho = cubit.create_vertex(  -7.9*km, -40.00*km, 0.0)
cubit.create_vertex(71.1*km, -20.00*km, 0.0)
p_slabtop_u = cubit.create_vertex(160.5*km, -7.50*km, 0.0)

cubit.cmd(f"create curve spline vertex {p_slabtop_west.id()} to {p_slabtop_u.id()} {p_topo_trench.id()}")
cubit.cmd(f"curve {cubit.get_last_id('curve')} name 'c_slabtop'")

# Bottom of slab (translate top of slab to the east)
#
# Better approach would be to move points normal to slab to preserve
# uniform thickness.
cubit.cmd(f"vertex {p_slabtop_west.id()} to {p_slabtop_coseismic.id()} copy move X {120.0*km}")
p_slabbot_west_id = p_slabtop_u.id() + 1

cubit.create_vertex(175.6*km, -40.0*km, 0.0)
p_slabbot_east = cubit.create_vertex(600.0*km, -40.0*km, 0.0)

cubit.cmd(f"create curve spline vertex {p_slabbot_west_id} to {p_slabbot_east.id()}")
cubit.cmd(f"curve {cubit.get_last_id('curve')} name 'c_slabbot'")

# Bottom edge of slab
cubit.cmd(f"create curve spline vertex {p_slabtop_west.id()} {p_slabbot_west_id}")
cubit.cmd(f"curve {cubit.get_last_id('curve')} name 'c_slab_end'")

# Top of mantle (uniform depth of 40 km)
p_moho_west = cubit.create_vertex(-600.0*km, -40.00*km, 0.0)

cubit.cmd(f"create curve spline vertex {p_moho_west.id()} {p_slabtop_moho.id()}")
cubit.cmd(f"curve {cubit.get_last_id('curve')} name 'c_conmoho'")

# Lateral edges and bottom boundary
p_bot_west = cubit.create_vertex(-600.0*km, -600.00*km, 0.0)
p_bot_east = cubit.create_vertex(+600.0*km, -600.00*km, 0.0)

cubit.cmd(f"create curve spline vertex {p_topo_west.id()} {p_moho_west.id()} {p_bot_west.id()}")
cubit.cmd(f"curve {cubit.get_last_id('curve')} name 'c_west'")

cubit.cmd(f"create curve spline vertex {p_bot_west.id()} {p_bot_east.id()}")
cubit.cmd(f"curve {cubit.get_last_id('curve')} name 'c_bot'")

cubit.cmd(f"create curve spline vertex {p_topo_east.id()} {p_slabbot_east.id()} {p_bot_east.id()}")
cubit.cmd(f"curve {cubit.get_last_id('curve')} name 'c_east'")


# Split curves to form bounding curves for each material
#
# Constructing the entire boundary curves as splines and then breaking
# them into pieces bounding the surfaces preserves continuity in slip.
cubit.cmd("split curve c_topo crossing curve c_slabtop")
cubit.cmd("split curve c_slabtop crossing curve c_conmoho")
cubit.cmd("split curve c_west crossing curve c_conmoho")
cubit.cmd("split curve c_east crossing curve c_slabbot")


# Create surfaces using bounding curves
# Continental crust
cubit.cmd("create surface curve c_topo c_west c_conmoho c_slabtop@A")
cubit.cmd(f"surface {cubit.get_last_id('surface')} name 'concrust'")

# Oceanic crust (slab)
cubit.cmd("create surface curve c_topo@A c_slabtop@A c_slabtop c_slab_end c_slabbot c_east")
cubit.cmd(f"surface {cubit.get_last_id('surface')} name 'oceancrust'")

# mantle
cubit.cmd("create surface curve c_west@A c_bot c_east@A c_slabbot c_slab_end c_slabtop c_conmoho")
cubit.cmd(f"surface {cubit.get_last_id('surface')} name 'mantle'")

# Imprint/merge
cubit.cmd("delete vertex all")
cubit.cmd("imprint all")
cubit.cmd("merge all")

# We must stitch the surfaces into a single volume in order to split
# the curves for the purpose of defining the discretization size along
# various portions of the curves.
cubit.cmd("stitch volume all")

# Split top of slab for fault surface
cubit.cmd(f"split curve c_slabtop distance {80.0*km} from end")

# Split topography/bathymetry to mimic same region as fault surface
# (used in setting discretization size)
cubit.cmd(f"split curve c_topo distance {190.0*km} from end")

# Split bottom of slab to mimic same region as fault surface
# (used in setting discretization size)
cubit.cmd(f"split curve c_slabbot distance {420.0*km} from end")
cubit.cmd(f"split curve c_slabbot distance {250.0*km} from end")


# -------------------------------------------------------------------------------------------------
# Generate the mesh
# -------------------------------------------------------------------------------------------------

cubit.cmd("surface all scheme trimesh")

# Reset sizes
cubit.cmd("curve all scheme default")
cubit.cmd("surface all sizing function none")

# Set size on faults (and bottom of slab corresponding to fault section)
cubit.cmd(f"curve c_slabtop@A c_slabtop size {DX_FAULT}")
cubit.cmd(f"curve c_topo size {DX_FAULT}")
cubit.cmd(f"curve c_slabbot@D size {DX_FAULT}")

# Use skeleton sizing to increase cell size away from fault.
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
cubit.cmd("block 1 surface concrust")
cubit.cmd("block 1 name 'concrust'")

cubit.cmd("block 2 surface oceancrust")
cubit.cmd("block 2 name 'oceancrust'")

cubit.cmd("block 3 surface mantle")
cubit.cmd("block 3 name 'mantle'")


# Nodesets

if False:
    # Fault (coseismic region only)
    cubit.cmd("group 'fault_coseismic' add node in c_slabtop@A")
    cubit.cmd("group 'fault_coseismic' add node in c_slabtop")
    cubit.cmd("nodeset 20 group fault_coseismic")
    cubit.cmd("nodeset 20 name 'fault_coseismic'")

    # Fault edge (coseismic region only)
    cubit.cmd("group 'fault_coseismic_edge' add node in vertex 65")
    cubit.cmd("nodeset 30 group fault_coseismic_edge")
    cubit.cmd("nodeset 30 name 'fault_coseismic_edge'")

    cubit.cmd("group 'fault_slabtop' add node in c_slabtop@A")
    cubit.cmd("group 'fault_slabtop' add node in c_slabtop")
    cubit.cmd("group 'fault_slabtop' add node in c_slabtop@D")
    cubit.cmd("nodeset 21 group fault_slabtop")
    cubit.cmd("nodeset 21 name 'fault_slabtop'")

    # Fault edge (slabtop)
    cubit.cmd("group 'fault_slabtop_edge' add node in vertex 33")
    cubit.cmd("nodeset 31 group fault_slabtop_edge")
    cubit.cmd("nodeset 31 name 'fault_slabtop_edge'")

    cubit.cmd("group 'fault_slabbot' add node in c_slabbot@B")
    cubit.cmd("group 'fault_slabbot' add node in c_slabbot@D")
    cubit.cmd("group 'fault_slabbot' add node in c_slabbot")
    cubit.cmd("nodeset 22 group fault_slabbot")
    cubit.cmd("nodeset 22 name 'fault_slabbot'")

    # Fault edge (slabbot)
    cubit.cmd("group 'fault_slabbot_edge' add node in vertex 24")
    cubit.cmd("nodeset 32 group fault_slabbot_edge")
    cubit.cmd("nodeset 32 name 'fault_slabbot_edge'")


    # Create nodeset for topography/bathymetry
    cubit.cmd("group 'groundsurf' add node in curve c_topo")
    cubit.cmd("group 'groundsurf' add node in curve c_topo@A")
    cubit.cmd("group 'groundsurf' add node in curve c_topo@B")
    cubit.cmd("nodeset 10 group groundsurf")
    cubit.cmd("nodeset 10 name 'groundsurf'")

    # Create nodesets for west boundary
    cubit.cmd("group 'bndry_west' add node in curve c_west")
    cubit.cmd("group 'bndry_west' add node in curve c_west@A")
    cubit.cmd("nodeset 11 group bndry_west")
    cubit.cmd("nodeset 11 name 'bndry_west'")

    # Create nodeset for east boundary
    # Crust
    cubit.cmd("group 'bndry_east_crust' add node in curve c_east")
    cubit.cmd("nodeset 12 group bndry_east_crust")
    cubit.cmd("nodeset 12 name 'bndry_east_crust'")

    # Mantle
    cubit.cmd("group 'bndry_east_mantle' add node in curve c_east@A")
    cubit.cmd("group 'bndry_east_mantle' remove node in group fault_slabbot")
    cubit.cmd("nodeset 13 group bndry_east_mantle")
    cubit.cmd("nodeset 13 name 'bndry_east_mantle'")

    # Create nodesets for bottom boundary
    cubit.cmd("group 'bndry_bot' add node in curve c_bot")
    cubit.cmd("nodeset 14 group bndry_bot")
    cubit.cmd("nodeset 14 name 'bndry_bot'")


# Sidesets
# Starting in PyLith v5, we will use sidesets instead of nodesets for BCs.

if True:
    # Fault (coseismic region only)
    cubit.cmd("group 'fault_coseismic' add c_slabtop@A")
    cubit.cmd("group 'fault_coseismic' add c_slabtop")
    cubit.cmd("sideset 20 group fault_coseismic")
    cubit.cmd("sideset 20 name 'fault_coseismic'")

    cubit.cmd("group 'fault_slabtop' add c_slabtop@A")
    cubit.cmd("group 'fault_slabtop' add c_slabtop")
    cubit.cmd("group 'fault_slabtop' add c_slabtop@D")
    cubit.cmd("sideset 21 group fault_slabtop")
    cubit.cmd("sideset 21 name 'fault_slabtop'")

    cubit.cmd("group 'fault_slabbot' add c_slabbot@B")
    cubit.cmd("group 'fault_slabbot' add c_slabbot@D")
    cubit.cmd("group 'fault_slabbot' add c_slabbot")
    cubit.cmd("sideset 22 group fault_slabbot")
    cubit.cmd("sideset 22 name 'fault_slabbot'")


    ## Create sideset for topography/bathymetry
    cubit.cmd("group 'groundsurf' add curve c_topo")
    cubit.cmd("group 'groundsurf' add curve c_topo@A")
    cubit.cmd("group 'groundsurf' add curve c_topo@B")
    cubit.cmd("sideset 10 group groundsurf")
    cubit.cmd("sideset 10 name 'groundsurf'")

    ## Create sideset for west boundary
    cubit.cmd("group 'bndry_west' add curve c_west")
    cubit.cmd("group 'bndry_west' add curve c_west@A")
    cubit.cmd("sideset 11 group bndry_west")
    cubit.cmd("sideset 11 name 'bndry_west'")

    # Create sideset for east boundary
    # Crust
    cubit.cmd("group 'bndry_east_crust' add curve c_east")
    cubit.cmd("sideset 12 group bndry_east_crust")
    cubit.cmd("sideset 12 name 'bndry_east_crust'")

    # Mantle
    cubit.cmd("group 'bndry_east_mantle' add curve c_east@A")
    cubit.cmd("group 'bndry_east_mantle' remove group fault_slabbot")
    cubit.cmd("sideset 13 group bndry_east_mantle")
    cubit.cmd("sideset 13 name 'bndry_east_mantle'")

    # Create sideset for bottom boundary
    cubit.cmd("group 'bndry_bot' add curve c_bot")
    cubit.cmd("sideset 14 group bndry_bot")
    cubit.cmd("sideset 14 name 'bndry_bot'")


    # Buried edges :TODO: Remove after using transform to create fault
    # Fault edge (coseismic region only)
    cubit.cmd("group 'fault_coseismic_edge' add node in vertex 65")
    cubit.cmd("nodeset 30 group fault_coseismic_edge")
    cubit.cmd("nodeset 30 name 'fault_coseismic_edge'")

    # Fault edge (slabtop)
    cubit.cmd("group 'fault_slabtop_edge' add node in vertex 33")
    cubit.cmd("nodeset 31 group fault_slabtop_edge")
    cubit.cmd("nodeset 31 name 'fault_slabtop_edge'")

    # Fault edge (slabbot)
    cubit.cmd("group 'fault_slabbot_edge' add node in vertex 24")
    cubit.cmd("nodeset 32 group fault_slabbot_edge")
    cubit.cmd("nodeset 32 name 'fault_slabbot_edge'")


# Write mesh as ExodusII file
cubit.cmd(f"export mesh 'mesh_{CELL}.exo' dimension 2 overwrite")
