#!python
# CUBIT script to generate simple Exodus files for testing ExodusII
# and MeshIOCubit.

import cubit

cubit.cmd("reset")

DX = 1.0 # mesh size

# Create geometry
v1 = cubit.create_vertex(-1.0,  0.0,  0.0)
v2 = cubit.create_vertex(+1.0,  0.0,  0.0)
v3 = cubit.create_vertex(+1.0, +1.0,  0.0)
v4 = cubit.create_vertex(-1.0, +1.0,  0.0)

c1 = cubit.create_curve(v1, v2)
c2 = cubit.create_curve(v2, v3)
c3 = cubit.create_curve(v3, v4)
c4 = cubit.create_curve(v4, v1)

s1 = cubit.create_surface([c1, c2, c3, c4])

# Mesh surface
cubit.cmd(f"surface {s1.id()} scheme trimesh")
cubit.cmd(f"surface {s1.id()} size {DX}")
cubit.cmd(f"mesh surface {s1.id()}")

# Create blocks
cubit.cmd("block 1 add tri 1 4")
cubit.cmd("block 1 name 'left'")

cubit.cmd("block 2 add tri 2 3")
cubit.cmd("block 2 name 'right'")

# Create nodesets
cubit.cmd(f"group 'vertices_xneg' add node in curve {c4.id()}")
cubit.cmd("nodeset 10 group vertices_xneg")
cubit.cmd("nodeset 10 name 'vertices_xneg'")

cubit.cmd(f"group 'vertices_yneg' add node in curve {c1.id()}")
cubit.cmd("nodeset 11 group vertices_yneg")
cubit.cmd("nodeset 11 name 'vertices_yneg'")

# create sidesets
cubit.cmd(f"group 'boundary_xneg' add curve {c4.id()}")
cubit.cmd("sideset 20 group boundary_xneg")
cubit.cmd("sideset 20 name 'boundary_xneg'")

cubit.cmd(f"group 'boundary_yneg' add curve {c1.id()}")
cubit.cmd("sideset 21 group boundary_yneg")
cubit.cmd("sideset 21 name 'boundary_yneg'")

version = cubit.get_version().split(".")[0]
filename = f"smalltri_v{version}.exo"
cubit.cmd(f"export mesh '{filename}' dimension 2 overwrite")
