#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2018 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/cornercases/3d/create_rotated_meshes.py
#
# @brief Create meshes rotated about (0,0,0)
#

import numpy
import math
from spatialdata.geocoords.CSCart import CSCart
from pylith.meshio.MeshIOAscii import MeshIOAscii

# Mesh filenames.
inputMeshes = ['onecell_hex.mesh', 'fivecells_tet.mesh']
outputMeshes = ['onecell_rotated_hex.mesh', 'fivecells_rotated_tet.mesh']
outFile = 'mesh_coords_rotated.txt'

# Just use mesh coordinates for now. I get an MPI error when I try to read the meshes.
meshCoords = numpy.array([[-4.0e+3, -4.0e+3,  0.0e+3,],
                          [-4.0e+3, -4.0e+3, -8.0e+3,],
                          [-4.0e+3, +4.0e+3, -8.0e+3,],
                          [-4.0e+3, +4.0e+3,  0.0e+3,],
                          [+4.0e+3, -4.0e+3,  0.0e+3,],
                          [+4.0e+3, -4.0e+3, -8.0e+3,],
                          [+4.0e+3, +4.0e+3, -8.0e+3,],
                          [+4.0e+3, +4.0e+3,  0.0e+3,]], dtype=numpy.float64)

# Rotation angles.
alpha = numpy.radians(30.0) # Rotation about x axis.
beta = numpy.radians(30.0) # Rotation about y axis.
# meshOrigin = numpy.array([0.0, 0.0, -4000.0], dtype=numpy.float64)

# Rotation matrix.
ca = math.cos(alpha)
sa = math.sin(alpha)
cb = math.cos(beta)
sb = math.sin(beta)
rotMatA = numpy.array([[ ca, -sa, 0.0],
                       [ sa,  ca, 0.0],
                       [0.0, 0.0, 1.0]], dtype=numpy.float64)
rotMatB = numpy.array([[ cb, 0.0,  sb],
                       [0.0, 1.0, 0.0],
                       [-sb, 0.0,  cb]], dtype=numpy.float64)
rotMat = numpy.dot(rotMatA, rotMatB)

# Loop over meshes.
"""
for meshNum in range(len(inputMeshes)):
    cs = CSCart()
    cs._configure()
    io = MeshIOAscii()
    inFile = inputMeshes[meshNum]
    outFile = outputMeshes[meshNum]
    io.inventory.filename = inFile
    io.inventory.coordsys = cs
    io._configure()

    mesh = io.read(debug=False)
"""
# meshCoords -= meshOrigin
coordsRot = numpy.dot(rotMat, meshCoords.transpose()).transpose()
numVerts = coordsRot.shape[0]
f = open(outFile, 'w')
fmt = '%d' + 3*'  %15.8f' + '\n'
for vertNum in range(numVerts):
    f.write(fmt % (vertNum, coordsRot[vertNum,0], coordsRot[vertNum,1], coordsRot[vertNum,2]))

f.close() 


# End of file

