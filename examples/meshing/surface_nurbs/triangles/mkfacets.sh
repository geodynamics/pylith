#! /bin/sh
# Shell script to create a triangulated surface in Facets format from a
# set of points. The Facets file can be read by Cubit.
# This script uses the GMT triangulate command, which invokes the triangle
# meshing package if you have installed it.

# Define filenames.
vertfile=../dem/ruapehu-nzmg-1km.txt
connectfile=ruapehu-nzmg-1km.connect
outfile=ruapehu-nzmg-1km.fac
vertidfile=ruapehu-nzmg-1km-id.vert
connectidfile=ruapehu-nzmg-1km-id.connect

# Triangulate the points to get connectivities.
triangulate $vertfile -V > $connectfile

# Count the number of vertices and cells to create a header.
numverts=$(wc -l < $vertfile)
numcells=$(wc -l < $connectfile)
header="$numverts $numcells"

# Create numbered versions of vertex and connectivity files (0-based).
nl -p -v0 $vertfile > $vertidfile
nl -p -v0 $connectfile > $connectidfile

# Cat files together to create Facets file.
echo $header > $outfile
cat $vertidfile $connectidfile >> $outfile
