#!python
# WARNING: This is not a standalone Python script. It must be run from within Cubit.
#
# Generate simple Exodus files for testing ExodusII and MeshIOCubit.
#
# pylith_convertmesh --reader=pylith.meshio.MeshIOCubit --reader.filename=tet_k.exo --writer=pylith.meshio.MeshIOAscii --writer.filename=tet_k.mesh


import cubit

# ----------------------------------------------------------------------
# Create block
# ----------------------------------------------------------------------

# Block is 40 km x 80 km x 40 km
# -20 km <= x <= 20 km
# -40 km <= y <= 40 km
# -40 km <= z <= 0 km
cubit.reset()

km = 1000.0
dx = 10*km

cubit.brick(40*km, 80*km, 40*km); v1 = cubit.get_last_id("volume")

# Translate block so the top is at z=0
cubit.cmd(f"volume {v1} move x 0 y 0 z -20*km")

# ----------------------------------------------------------------------
# Create interface surfaces
# ----------------------------------------------------------------------
cubit.cmd("create planar surface with plane xplane offset 0"); fault_surf = cubit.get_last_id("surface")
v_fault = cubit.get_last_id("volume")
cubit.cmd(f"surface {fault_surf} name 'fault_surface'")

# ----------------------------------------------------------------------
# Divide volumes using interface surfaces
# ----------------------------------------------------------------------
cubit.cmd(f"webcut volume {v1} with plane surface {fault_surf}"); v2 = cubit.get_last_id("volume")
cubit.cmd(f"volume {v1} name 'v_xpos'"); v_xpos = v1
cubit.cmd(f"volume {v2} name 'v_xneg'"); v_xneg = v2

# ----------------------------------------------------------------------
# Chop x and y faults into pieces
# ----------------------------------------------------------------------
# fault_surface_x
cubit.brick(2*km, 50*km, 40*km); fault_chopper = cubit.get_last_id("volume")
cubit.cmd(f"chop volume {v_fault} with volume {fault_chopper}")

# ----------------------------------------------------------------------
# Imprint all volumes, merging surfaces
# ----------------------------------------------------------------------
cubit.cmd("imprint all with volume all")
cubit.cmd("merge all")
cubit.cmd("delete body 5 6")

# ----------------------------------------------------------------------
# Set discretization size
# ----------------------------------------------------------------------
cubit.cmd(f"volume all size {dx}")

# ----------------------------------------------------------------------
# Generate the mesh
# ----------------------------------------------------------------------
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

cubit.cmd("volume all smooth scheme condition number beta 1.7 cpu 10")
cubit.cmd("smooth volume all")

# ----------------------------------------------------------------------
# Create blocks for materials
# ----------------------------------------------------------------------
cubit.cmd(f"block 10 volume {v_xpos} {v_xneg}")
cubit.cmd("block 10 name 'elastic'")

# ----------------------------------------------------------------------
# Create nodesets for faults
# ----------------------------------------------------------------------
cubit.cmd(f"group 'fault' add node in surface fault_surface")
cubit.cmd("nodeset 10 group fault")
cubit.cmd("nodeset 10 name 'fault'")

cubit.cmd("group 'fault_edge' add node in curve 45")
cubit.cmd("group 'fault_edge' add node in curve 46")
cubit.cmd("group 'fault_edge' add node in curve 47")
cubit.cmd("nodeset 11 group fault_edge")
cubit.cmd("nodeset 11 name 'fault_edge'")

# ----------------------------------------------------------------------
# Create nodeset for +z face
# ----------------------------------------------------------------------
cubit.cmd("group 'vertices_zpos' add node in surface 9")
cubit.cmd("group 'vertices_zpos' add node in surface 16")
cubit.cmd("nodeset 20 group vertices_zpos")
cubit.cmd("nodeset 20 name 'vertices_zpos'")

# ----------------------------------------------------------------------
# Create sideset for fault
# ----------------------------------------------------------------------
cubit.cmd(f"group 'fault_faces' add surface fault_surface")
cubit.cmd("sideset 110 group fault_faces")
cubit.cmd("sideset 110 name 'fault_faces'")

# ----------------------------------------------------------------------
# Create nodeset for +z face
# ----------------------------------------------------------------------
cubit.cmd("group 'boundary_zpos' add surface 9")
cubit.cmd("group 'boundary_zpos' add surface 16")
cubit.cmd("sideset 120 group boundary_zpos")
cubit.cmd("sideset 120 name 'boundary_zpos'")

# ----------------------------------------------------------------------
# Export exodus file
# ----------------------------------------------------------------------
cubit.cmd("export mesh 'tet_k.exo' dimension 3 overwrite")
