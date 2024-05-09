#!python
# CUBIT script to generate simple Exodus files for testing ExodusII
# and MeshIOCubit.

import cubit

cubit.cmd("reset")

DX = 1.0 # mesh size

# Create geometry
v1 = cubit.brick(2.0, 1.0, 1.0)
s_xneg = 4
s_zneg = 2

# Mesh surface
cubit.cmd(f"volume {v1.id()} scheme map")
cubit.cmd(f"volume {v1.id()} size {DX}")
cubit.cmd(f"mesh volume {v1.id()}")

# Create blocks
cubit.cmd(f"block 1 add volume {v1.id()}")
cubit.cmd("block 1 name 'domain'")

# Create nodesets
cubit.cmd(f"group 'vertices_xneg' add node in surface {s_xneg}")
cubit.cmd("nodeset 10 group vertices_xneg")
cubit.cmd("nodeset 10 name 'vertices_xneg'")

cubit.cmd(f"group 'vertices_zneg' add node in surface {s_zneg}")
cubit.cmd("nodeset 11 group vertices_zneg")
cubit.cmd("nodeset 11 name 'vertices_zneg'")

# create sidesets
cubit.cmd(f"group 'boundary_xneg' add surface {s_xneg}")
cubit.cmd("sideset 20 group boundary_xneg")
cubit.cmd("sideset 20 name 'boundary_xneg'")

cubit.cmd(f"group 'boundary_zneg' add surface {s_zneg}")
cubit.cmd("sideset 21 group boundary_zneg")
cubit.cmd("sideset 21 name 'boundary_zneg'")

version = cubit.get_version().split(".")[0]
filename = f"smallhex_v{version}.exo"
cubit.cmd(f"export mesh '{filename}' dimension 3 overwrite")
