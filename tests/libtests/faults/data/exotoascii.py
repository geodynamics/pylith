#!/usr/bin/env nemesis

import sys

if len(sys.argv) != 2:
    raise ValueError("usage: exotoascii.py EXODUS_MESH")
filenameEXO = sys.argv[1]
filenameASCII = filenameEXO.replace(".exo", ".mesh")

from pylith.meshio.MeshIOAscii import MeshIOAscii
from pylith.meshio.MeshIOCubit import MeshIOCubit
import pylith.utils.petsc as petsc
from spatialdata.geocoords.CSCart import CSCart

petsc.initialize(sys.argv)

cs = CSCart()
cs._configure()
cubit = MeshIOCubit()
cubit.inventory.filename = filenameEXO
cubit.inventory.coordsys = cs
cubit._configure()
cubit.preinitialize()

mesh = cubit.read(check=True)

ascii = MeshIOAscii()
ascii.inventory.filename = filenameASCII
ascii._configure()
ascii.preinitialize()
ascii.write(mesh)

del cubit
del ascii
del mesh

petsc.finalize()

# End of file
