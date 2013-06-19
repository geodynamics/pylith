#!/usr/bin/env python
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# PREREQUISITES: numpy

# ======================================================================
import numpy
# import pdb
# pdb.set_trace()

# Define parameters.
splayFile = "splay_points2.txt"
numProfs = 5
pointsPerProf = 4
journalFile = "splay_skinsurf.jou"
acisFile = "splay_surf.sat"

# Journal file formatting, etc.
separator = "# ----------------------------------------------------------\n"
journalBeg = \
           "reset\n" + \
           "# This is a simple Cubit journal file to create an ACIS\n" + \
           "# NURBS surface from a set of profiles.\n" + \
           "#\n" + separator
lineBeg = "create curve spline"
splineFmt = " location %15.11e %15.11e %15.11e"
skinCmd = "create surface skin curve all\n"
delCmd = "delete curve all\n"
expCmd = "export Acis '" + acisFile + "'\n"


# Read coordinates and reshape them.
splayCoords = numpy.loadtxt(splayFile, dtype=numpy.float64).reshape(
    numProfs, pointsPerProf, 3)

j = open(journalFile, 'w')
j.write(journalBeg)

# Loop over profiles.
for profile in range(numProfs):
    points = splayCoords[profile,:,:]
    j.write(lineBeg)
    for pointNum in range(pointsPerProf):
        point = points[pointNum,:]
        j.write(splineFmt % (point[0], point[1], point[2]))

    j.write("\n")

# Create skin surface.
j.write(skinCmd)

# Delete spline curves and export Acis file.
j.write(delCmd)
j.write(expCmd)
j.close()

# End of file
