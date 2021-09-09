#!/usr/bin/env python
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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/manual/2d/powerlaw/readnormals.py

## @brief Read debugging output from cylinder test problem.

import math
import numpy

# ----------------------------------------------------------------------
# Filenames.
prefixes = ['cylinder_pres_elastic_norefstate_b1_q1_quad', 'cylinder_pres_elastic_norefstate_b1_q1_tri',
            'cylinder_pres_elastic_norefstate_b2_q2_quad', 'cylinder_pres_elastic_norefstate_b2_q2_tri',
            'cylinder_pres_elastic_norefstate_b3_q3_quad', 'cylinder_pres_elastic_norefstate_b3_q3_tri',
            'cylinder_pres_elastic_norefstate_b4_q4_quad', 'cylinder_pres_elastic_norefstate_b4_q4_tri',
            'cylinder_pres_elastic_norefstate_b1_q1_fine_quad', 'cylinder_pres_elastic_norefstate_b1_q1_fine_tri',
            'cylinder_pres_elastic_norefstate_b2_q2_fine_quad', 'cylinder_pres_elastic_norefstate_b2_q2_fine_tri',
            'cylinder_pres_elastic_norefstate_b3_q3_fine_quad', 'cylinder_pres_elastic_norefstate_b3_q3_fine_tri',
            'cylinder_pres_elastic_norefstate_b4_q4_fine_quad', 'cylinder_pres_elastic_norefstate_b4_q4_fine_tri']
string1 = '0 TS'
string2 = '0 SNES'
lengthScale = 1000.0
numLinesPerPoint = 6

# ----------------------------------------------------------------------
def getPoints(filePrefix):
    """
    Get info for a set of points.
    """
    inFile = filePrefix + '.log'
    f = open(inFile, 'r')
    lines = f.readlines()
    string1Lines = [line for line in lines if string1 in line]
    startLine = lines.index(string1Lines[0]) + 1
    string2Lines = [line for line in lines if string2 in line]
    endLine = lines.index(string2Lines[0])
    linesUse = lines[startLine:endLine]
    numLines = len(linesUse)
    numPoints = numLines//numLinesPerPoint
    x = numpy.zeros((numPoints, 3), dtype=numpy.float64)
    feNorm = numpy.zeros((numPoints, 3), dtype=numpy.float64)
    anlNorm = numpy.zeros((numPoints, 3), dtype=numpy.float64)
    feg0 = numpy.zeros((numPoints, 3), dtype=numpy.float64)
    anlg0 = numpy.zeros((numPoints, 3), dtype=numpy.float64)
    numRead = 0

    for line in linesUse:
        ind = line.find(':')
        if (ind != -1):
            val = line[:ind]
            lineStrip = line[ind+1:].split()
            if (val == 'x'):
                x[numRead,:2] = [float(i) for i in lineStrip]
            elif (val == 'FE normal'):
                feNorm[numRead,:2] = [float(i) for i in lineStrip]
            elif (val == 'Anl normal'):
                anlNorm[numRead,:2] = [float(i) for i in lineStrip]
            elif (val == 'FE g0'):
                feg0[numRead,:2] = [float(i) for i in lineStrip]
            elif (val == 'Anl g0'):
                anlg0[numRead,:2] = [float(i) for i in lineStrip]
                numRead += 1

    f.close()

    x *= lengthScale

    writeVtk(filePrefix, x, feNorm, anlNorm, feg0, anlg0)

    return

def writeVtk(filePrefix, x, feNorm, anlNorm, feg0, anlg0):
    """
    Write VTK file for a set of points.
    """
    numPoints = x.shape[0]
    fileName = filePrefix + '.vtk'
    vtkHead = '# vtk DataFile Version 2.0\n' + \
        'Kernel information on Neumann boundary\n' + \
        'ASCII\n' + \
        'DATASET POLYDATA\n' + \
        'POINTS %d double\n' % numPoints
    v = open(fileName, 'w')
    v.write(vtkHead)
    numpy.savetxt(v, x)
    pointHead = 'VERTICES %d %d\n' % (numPoints, 2*numPoints)
    verts = numpy.arange(numPoints, dtype=numpy.int64)
    vertsPerCell = numpy.ones_like(verts)
    outVerts = numpy.column_stack((vertsPerCell, verts))
    v.write(pointHead)
    numpy.savetxt(v, outVerts, fmt='%d')

    vtkHead2 = 'POINT_DATA %d\n' % numPoints
    vtkHead2a = 'VECTORS FE_normal double\n'
    v.write(vtkHead2)
    v.write(vtkHead2a)
    numpy.savetxt(v, feNorm)

    vtkHead3 = 'VECTORS ANL_normal double\n'
    v.write(vtkHead3)
    numpy.savetxt(v, anlNorm)

    vtkHead4 = 'VECTORS FE_g0 double\n'
    v.write(vtkHead4)
    numpy.savetxt(v, feg0)

    vtkHead5 = 'VECTORS ANL_g0 double\n'
    v.write(vtkHead5)
    numpy.savetxt(v, anlg0)

    v.close()

    return
# ----------------------------------------------------------------------

for filePrefix in prefixes:
    getPoints(filePrefix)
    
# End of file 
