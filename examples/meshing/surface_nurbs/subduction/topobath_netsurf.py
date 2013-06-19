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
demFile = "topobath_points.txt"
numProfs = 4
pointsPerProf = 9
journalFile = "topobath_netsurf.jou"
acisFile = "topobath_surf.sat"

# Journal file formatting, etc.
separator = "# ----------------------------------------------------------\n"
journalBeg = \
           "reset\n" + \
           "# This is a simple Cubit journal file to create an ACIS\n" + \
           "# NURBS surface from a set of intersecting lines.\n" + \
           "#\n" + separator
lineBeg = "create curve spline"
splineFmt = " location %15.11e %15.11e %15.11e"
netFmt = "create surface net u curve %d to %d v curve %d to %d\n"
delCmd = "delete curve all\n"
expCmd = "export Acis '" + acisFile + "'\n"


# Read coordinates and reshape them.
demCoords = numpy.loadtxt(demFile, dtype=numpy.float64).reshape(
    numProfs, pointsPerProf, 3)

j = open(journalFile, 'w')
j.write(journalBeg)

# Loop over profiles (u-lines).
for profile in range(numProfs):
    points = demCoords[profile,:,:]
    j.write(lineBeg)
    for pointNum in range(pointsPerProf):
        point = points[pointNum,:]
        j.write(splineFmt % (point[0], point[1], point[2]))

    j.write("\n")

# Loop over contours (v-lines).
for contour in range(pointsPerProf):
    points = demCoords[:,contour,:]
    j.write(lineBeg)
    for pointNum in range(numProfs):
        point = points[pointNum,:]
        j.write(splineFmt % (point[0], point[1], point[2]))

    j.write("\n")

# Create net surface.
uline1 = 1
uline2 = numProfs
vline1 = numProfs + 1
vline2 = vline1 + pointsPerProf - 1
j.write(netFmt % (uline1, uline2, vline1, vline2))

# Delete spline curves and export Acis file.
j.write(delCmd)
j.write(expCmd)
j.close()

# End of file
