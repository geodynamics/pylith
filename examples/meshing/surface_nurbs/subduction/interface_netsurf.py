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
intFile = "interface_points.txt"
numConts = 5
pointsPerCont = 5
journalFile = "interface_netsurf.jou"
acisFile = "interface_surf.sat"

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
intCoords = numpy.loadtxt(intFile, dtype=numpy.float64).reshape(
    numConts, pointsPerCont, 3)

j = open(journalFile, 'w')
j.write(journalBeg)

# Loop over contours (u-lines).
for contour in range(numConts):
    points = intCoords[contour,:,:]
    j.write(lineBeg)
    for pointNum in range(pointsPerCont):
        point = points[pointNum,:]
        j.write(splineFmt % (point[0], point[1], point[2]))

    j.write("\n")

# Loop over profiles (v-lines).
for profile in range(pointsPerCont):
    points = intCoords[:,profile,:]
    j.write(lineBeg)
    for pointNum in range(numConts):
        point = points[pointNum,:]
        j.write(splineFmt % (point[0], point[1], point[2]))

    j.write("\n")

# Create net surface.
uline1 = 1
uline2 = numConts
vline1 = numConts + 1
vline2 = vline1 + pointsPerCont - 1
j.write(netFmt % (uline1, uline2, vline1, vline2))

# Delete spline curves and export Acis file.
j.write(delCmd)
j.write(expCmd)
j.close()

# End of file
