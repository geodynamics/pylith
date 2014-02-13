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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# PREREQUISITES: numpy

# ======================================================================
import numpy

# Define parameters.
splayFile = "splay_points.txt"
numProfiles = 5
pointsPerProfile = 4
journalFile = "splay_skinsurf.jou"
acisFile = "splay_surf.sat"

# Journal file formatting, etc.
separator = "# ----------------------------------------------------------\n"
journalBeg = \
           "# CUBIT journal file generated by splay_skinsurf.py.\n" + \
           "#\n" + \
           "# Create an ACIS NURBS surface from a set of profiles.\n" + \
           "#\n" + separator + \
           "reset\n"
lineBeg = "create curve spline"
splineFmt = " location %10.2e %10.2e %10.2e"
skinCmd = "create surface skin curve all\n"
delCmd = "delete curve all\n"
expCmd = "export Acis '" + acisFile + "'\n"


# Read coordinates and reshape them.
splayCoords = numpy.loadtxt(splayFile, dtype=numpy.float64).reshape(numProfiles, pointsPerProfile, 3)

j = open(journalFile, 'w')
j.write(journalBeg)

# Loop over profiles.
for profile in range(numProfiles):
    points = splayCoords[profile,:,:]
    j.write(lineBeg)
    for pointNum in range(pointsPerProfile):
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
