#!/usr/bin/env nemesis

import sys

if len(sys.argv) != 2:
    raise ValueError("usage: exotoascii.py EXODUS_MESH")
filenameEXO = sys.argv[1]
filenameASCII = filenameEXO.replace(".exo", ".mesh")

from pylith.meshio.MeshIOAscii import MeshIOAscii
from pylith.meshio.MeshIOCubit import MeshIOCubit
import pylith.utils.petsc as petsc

petsc.initialize(sys.argv)

cubit = MeshIOCubit()
cubit.inventory.filename = filenameEXO
cubit._configure()

mesh = cubit.read(debug=False, interpolate=True)

ascii = MeshIOAscii()
ascii.inventory.filename = filenameASCII
ascii._configure()
ascii.write(mesh)

del cubit
del ascii
del mesh

petsc.finalize()

# End of file
