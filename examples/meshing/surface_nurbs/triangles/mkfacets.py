#!/usr/bin/env python
# Python script to create a triangulated surface in Facets format from a
# set of points. The Facets file can be read by Cubit.
# This script uses the GMT triangulate command, which invokes the triangle
# meshing package if you have installed it.

# Import necessary packages.
import os
import os
# import pdb
# pdb.set_trace()

# Define filenames.
vertexFile = "../dem/ruapehu-nzmg-1km.txt"
connectFile = "ruapehu-nzmg-1km.connect"
facetsOut = "ruapehu-nzmg-1km.fac"

# Triangulate the points to get connectivities.
os.system('triangulate '+vertexFile+' > '+connectFile)

# Read connectivities and vertices.
v = open(vertexFile, 'r')
vertices = v.readlines()
v.close()
numVertices = len(vertices)
c = open(connectFile, 'r')
connect = c.readlines()
c.close()
numTriangles = len(connect)

# Create factes file header and write it to facets file.
header = "%d  %d\n" % (numVertices, numTriangles)
f = open(facetsOut, 'w')
f.write(header)

# Write numbered vertices to facets file.
vertFormat = "%d  %s"
for vertex in range(numVertices):
  outLine = vertFormat % (vertex, vertices[vertex])
  f.write(outLine)
  
# Write numbered connectivities to facets file.
connectFormat = "%d  %s"
for triangle in range(numTriangles):
  outLine = connectFormat % (triangle, connect[triangle])
  f.write(outLine)
f.close()
